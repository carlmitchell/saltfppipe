from sys import exit as crash
import numpy as np
import matplotlib.pyplot as plt
from photutils import daofind
from saltfppipe.fp_image_class import FPImage

class Verify_Center_Plot:
    def __init__(self,fn,data,flag,xcen,ycen,objx,objy,ghox,ghoy):
        self.key = None
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #Display the image data and set up the plot object
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
        #Connect the dots for the object/ghost matches
        for i in range(len(objx)):
            plt.plot([objx[i],ghox[i]],[objy[i],ghoy[i]],color="blue")
        #Mark the center
        plt.scatter([xcen],[ycen],color="red",s=30)
        plt.xlim(0,data.shape[1])
        plt.ylim(0,data.shape[0])
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",
                                                       self.keypress)
        plt.tight_layout()
        plt.show()
        
    def keypress(self,event):
        if event.key in ["w","a","s","d","q"]:
            self.key = event.key
            plt.close()

def verify_center(fnlist, objxs, objys, ghoxs, ghoys, xcens, ycens):
    """Displays images and their likely ghost pairs and allows the user to
    call a crash if the center locations aren't good.
    
    """
    
    #Change pyplot shortcut keys
    oldsavekey = plt.rcParams["keymap.save"]
    plt.rcParams["keymap.save"] = ""
    
    #Make empty lists to load things for plotting
    flagged = np.zeros(len(fnlist),dtype="bool")
    
    #Open each image, make lists of the objects and ghosts in each image
    images = [FPImage(fn) for fn in fnlist]
    axcens = [image.axcen for image in images]
    aycens = [image.aycen for image in images]
    arads = [image.arad for image in images]
    datas = [image.inty[aycen-arad:aycen+arad, axcen-arad:axcen+arad]
             for (image,axcen,aycen,arad) in zip(images,axcens,aycens,arads)]
    
    #Adjust object ghost, and center coordinates because of the aperture mask
    for i in range(len(fnlist)):
        objxs[i] = np.array(objxs[i])-axcens[i]+arads[i]
        objys[i] = np.array(objys[i])-aycens[i]+arads[i]
        ghoxs[i] = np.array(ghoxs[i])-axcens[i]+arads[i]
        ghoys[i] = np.array(ghoys[i])-aycens[i]+arads[i]
        xcens[i] = xcens[i]-axcens[i]+arads[i]
        ycens[i] = ycens[i]-aycens[i]+arads[i]
    
    #Interactive plotting routine thing
    i = 0
    while True:
        plot = Verify_Center_Plot(fnlist[i],datas[i],
                                  flagged[i],xcens[i],ycens[i],
                                  objxs[i],objys[i],ghoxs[i],ghoys[i])
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
    for image in images: image.close()
    
    #Restore old save key
    plt.rcParams["keymap.save"] = oldsavekey
    
    #If any images are flagged, return False
    if np.sum(flagged)>0:
        print repr(np.int(np.sum(flagged)))+" centers rejected!"
        return False
    else:
        print "All centers approved!"
        return True
    
def find_ghost_centers(fnlist, tolerance=3, thresh=4.5, guess=None):
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
    guess -> Optional. If you already have an idea of where the center should
             be, you'd put it here. Should be a 2-long iterable, with guess[0]
             being the X center guess, and guess[1] being the y center guess.

    Outputs:
    xcenlist -> List of image center X coordinates
    ycenlist -> List of image center Y coordinates
    
    """
    
    #Get image FWHMs
    fwhm = np.empty(len(fnlist))
    firstimage = FPImage(fnlist[0])
    toggle = firstimage.fwhm
    axcen = firstimage.axcen
    aycen = firstimage.aycen
    arad = firstimage.arad
    firstimage.close()
    if axcen == None:
        print "Error! Images have not yet been aperture-masked! Do this first!"
        crash()
    if toggle == None:
        print "Warning: FWHMs have not been measured!"
        print "Assuming 5 pixel FWHM for all images."
        fwhm = 5.*np.ones(len(fnlist))
    else:
        for i in range(len(fnlist)):
            image = FPImage(fnlist[i])
            fwhm[i] = image.fwhm
            image.close()
            
    #Get sky background levels
    skyavg = np.empty(len(fnlist))
    skysig = np.empty(len(fnlist))
    for i in range(len(fnlist)):
        image = FPImage(fnlist[i])
        skyavg[i], skysig[i], _skyvar = image.skybackground()
        image.close()
        
    #Identify the stars in each image
    xlists = []
    ylists = []
    maglists = []
    print "Identifying stars and ghosts in each image..."
    for i in range(len(fnlist)):
        xlists.append([])
        ylists.append([])
        maglists.append([])
        image = FPImage(fnlist[i])
        axcen = image.axcen
        aycen = image.aycen
        arad = image.arad
        sources = daofind(image.inty-skyavg[i],
                          fwhm=fwhm[i],
                          threshold=thresh*skysig[i]).as_array()
        for j in range(len(sources)):
            #Masks for center and edge of image
            cenmask = ((sources[j][1]-axcen)**2 +
                       (sources[j][2]-aycen)**2 > (0.05*arad)**2)
            edgemask = ((sources[j][1]-axcen)**2 +
                        (sources[j][2]-aycen)**2 < (0.95*arad)**2)
            if np.logical_and(cenmask, edgemask):
                xlists[i].append(sources[j][1])
                ylists[i].append(sources[j][2])
                maglists[i].append(sources[j][-1])
        image.close()
    
    if guess is None:
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
                dist2array = ( (xcenarray-xcenarray[j])**2 +
                               (ycenarray-ycenarray[j])**2 )
                sub1array[j] = np.sum(dist2array<tolerance**2)
            #Determine the locations of the "best" centers.
            bestcenloc = np.where(sub1array==max(sub1array))[0]
            #Now cross check JUST the best ones against each other
            sub1array = np.zeros(len(bestcenloc))
            xcenarray = xcenarray[bestcenloc]
            ycenarray = ycenarray[bestcenloc]
            for j in range(len(bestcenloc)):
                dist2array = ( (xcenarray-xcenarray[j])**2 +
                               (ycenarray-ycenarray[j])**2 )
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
        bestxcen = np.average(xcen[bestcenloc])
        bestycen = np.average(ycen[bestcenloc])
    
    else:
        #Forced guess:
        bestxcen, bestycen = guess[0], guess[1]
    
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
            dist2list = ( (xlists[i] - refxlist[j])**2 + 
                          (ylists[i] - refylist[j])**2 )
            matchlist = dist2list < tolerance**2
            if np.sum(matchlist) >= 1:
                #We found a match! Now we need to know where the match is
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
        image = FPImage(fnlist[i],update=True)
        image.xcen = xcenlist[i]
        image.ycen = ycenlist[i]
        image.close()
        
    #Manually verify ghost centers
    while True:
        yn = raw_input("Manually verify ghost centers? (Recommended) (y/n) ")
        if "n" in yn or "N" in yn:
            break
        elif "y" in yn or "Y" in yn:
            goodtog = verify_center(fnlist,
                                    objxs,objys,
                                    ghoxs,ghoys,
                                    xcenlist,ycenlist)
            if goodtog: break
            else:
                print "Centers not approved!"
                crash()
        
    return xcenlist, ycenlist
