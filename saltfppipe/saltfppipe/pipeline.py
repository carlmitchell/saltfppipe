import sys
from os.path import isdir, isfile, join, split
from os import mkdir, listdir, remove
from shutil import rmtree, copyfile, move
import numpy as np
import matplotlib.pyplot as plt
from astropy.io.fits import open as openfits
from scipy.optimize import minimize
from saltfppipe.aperture_mask import aperture_mask
from saltfppipe.combine_flat import combine_flat
from saltfppipe.flatten import flatten
from saltfppipe.find_ghost_centers import find_ghost_centers
#from lacosmicx import lacosmicx
from saltfppipe.deghost import deghost
from saltfppipe.align_norm import align_norm
from saltfppipe.fit_wave_soln import fit_wave_soln
from astropy.io.fits import PrimaryHDU
from saltfppipe.gauss_fit import GaussFit
from saltfppipe.sub_sky_rings import sub_sky_rings
from saltfppipe.make_final_image import make_final_image
from saltfppipe.solar_velocity_shift import solar_velocity_shift
from saltfppipe.fit_velmap import fit_velmap_ha_n2_mode
from pyraf import iraf
from zsalt.imred import imred
from saltfppipe.make_clean_map import make_clean_map

class Verify_Images_Plot:
    def __init__(self,data,obj_type,flagged,num,total):
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
        plt.title("Image #"+str(num)+"/"+str(total)+", Press 'q' to quit.\n"+
                  "Object = "+obj_type+", Flagged = "+str(flagged)+"\n"+
                  "Use 'a'/'d' Keys to Navigate \n"+
                  "'w' and 's' Keys to Flag Bad Images/Headers")
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
        
    def keypress(self, event):
        if event.key in ["a", "s", "w", "d", "q"]:
            self.key = event.key
            plt.close()
        return

class Review_Images_Plot:
    def __init__(self,data,obj_type,filename):
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
        plt.title("File: "+filename+", Object: "+obj_type+"\n"+
                  "Press 'h' to correct header if object type wrong.\n"+
                  "Press 'd' to delete this image.\n"+
                  "Press 'a' to approve this image.")
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
        
    def keypress(self, event):
        if event.key in ["h","d","a"]:
            self.key = event.key
            plt.close()

class Click_Near_Aperture_Plot:
    def __init__(self,data,xs,ys,axcen,aycen,arad):
        #Initialize xcoo and ycoo
        self.xcoo = None
        self.ycoo = None
        #Create a figure
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #Plot the image
        self.plot = ax1.imshow(data,cmap="Greys",aspect="equal",vmin=np.percentile(data,1),vmax=np.percentile(data,99),origin="lower",interpolation="nearest")
        #Plot the points on the boundary
        ax1.scatter(xs,ys,color="red")
        #Plot the best fitting circle center
        ax1.scatter(axcen,aycen,color="blue")
        #Plot the best fitting circle center
        xcirc = axcen+arad*np.cos(np.linspace(0,2*np.pi,1000))
        ycirc = aycen+arad*np.sin(np.linspace(0,2*np.pi,1000))
        ax1.plot(xcirc,ycirc,color="blue")
        #Fix the extent of the plot
        plt.xlim(0,data.shape[1])
        plt.ylim(0,data.shape[0])
        
        plt.title("Left click to zoom in on a point near the aperture boundary.\n"+
                  "Right click to quit.")
        self.cid = self.plot.figure.canvas.mpl_connect("button_press_event",self.click)
        plt.tight_layout()
        plt.show()
    
    def click(self,event):
        if event.button == 1:
            if event.inaxes != self.plot.axes: return
            self.xcoo = event.xdata
            self.ycoo = event.ydata
            plt.close()
            return
        else:
            plt.close()
            return

class Zoom_Aperture_Plot:
    def __init__(self,data):
        #Initialize xcoo and ycoo
        self.xcoo = None
        self.ycoo = None
        #Create a figure
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #Plot the image
        self.plot = ax1.imshow(data,cmap="Greys",aspect="equal",vmin=np.percentile(data,1),vmax=np.percentile(data,99),origin="lower",interpolation="nearest")
        plt.title("Left click to mark a boundary location.\n"+
                  "Right click to go back.")
        plt.tight_layout()
        self.cid = self.plot.figure.canvas.mpl_connect("button_press_event",self.click)
        plt.show()
    
    def click(self,event):
        if event.button == 1:
            if event.inaxes != self.plot.axes: return
            self.xcoo = event.xdata
            self.ycoo = event.ydata
            plt.close()
            return
        else:
            plt.close()
            return

class Mark_Stars_Plot:
    def __init__(self,data,xs,ys):
        #Initialize x and y coordinates
        self.xcoo = None
        self.ycoo = None
        self.key = None
        #Create a figure
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #Plot the image
        self.plot = ax1.imshow(data,cmap="Greys",aspect="equal",vmin=np.percentile(data,1),vmax=np.percentile(data,99),origin="lower",interpolation="nearest")
        #Plot the previously marked stars
        ax1.scatter(xs,ys,color="red")
        #Fix the extent of the plot
        self.extent = data.shape
        plt.xlim(0,self.extent[1])
        plt.ylim(0,self.extent[0])
        plt.title("Press Z to zoom in on a region.\n"+
                  "Press D to delete a marked star.\n"+
                  "Press Q to quit.")
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
    
    def keypress(self,event):
        self.key = event.key
        if self.key == "q":
            plt.close()
            return
        if self.key in ["d","z"] and event.inaxes == self.plot.axes and event.xdata>24 and event.xdata<self.extent[1]-25 and event.ydata>24 and event.ydata<self.extent[0]-25:
            self.xcoo = event.xdata
            self.ycoo = event.ydata
            plt.close()
            return
    
class Zoom_Mark_Stars_Plot:
    def __init__(self,data):
        #Initialize x and y coordinates
        self.xcoo = None
        self.ycoo = None
        #Create figure
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #Plot the image
        self.plot = ax1.imshow(data,cmap="Greys",aspect="equal",vmin=np.percentile(data,1),vmax=np.percentile(data,99),origin="lower",interpolation="nearest")
        #Fix the extent of the plot
        plt.xlim(0,data.shape[1])
        plt.ylim(0,data.shape[0])
        plt.title("Press Z to zoom back out.\n"+
                  "Press W to mark a star.")
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.show()
    
    def keypress(self,event):
        if event.key == "z":
            plt.close()
            return
        if event.key == "w" and event.inaxes == self.plot.axes:
            self.xcoo = event.xdata
            self.ycoo = event.ydata
            plt.close()
            return

class FWHMs_Plot:
    def __init__(self,fwhms,oldfwhms,oldfwhmserr):
        #Initialize button
        self.button = None
        #Create figure
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #Plot the fwhm data
        self.plot = ax1.scatter(np.arange(len(fwhms)),fwhms,color="blue",label="New FWHM")
        if not (oldfwhms is None) and not (oldfwhmserr is None):
            ax1.errorbar(np.arange(len(fwhms)),oldfwhms,yerr=oldfwhmserr,color="red",label="Previous FWHMs",fmt="o")
        if not (oldfwhms is None) and oldfwhmserr is None:
            ax1.scatter(np.arange(len(fwhms)),oldfwhms,color="red",label="Previous FWHM")
        plt.title("Left click to accept new star's measurements\n"+
                  "Right click to reject it.")
        plt.legend(loc="best")
        plt.xlim(-1,len(fwhms)+1)
        plt.xlabel("Image number in track")
        plt.ylabel("FWHM [pix]")
        self.cid = self.plot.figure.canvas.mpl_connect("button_press_event",self.click)
        plt.show()
    
    def click(self,event):
        if event.button == 1: self.button = 1
        else: self.button = None
        plt.close()
        return
    
def create_directory(name):
    """Procedure to create a new directory. If the directory already
    exists, asks user to overwrite or not. If the old directory is not
    overwritten, returns False. Otherwise returns True.

    Inputs:
    name -> String containing the name of the directory to be created
    
    Outputs:
    Boolean toggle, whether or not the directory was created

    """
    if isdir(name):
        while True:
            yn = raw_input("Directory '"+name+"' already exists. Overwrite it? (y/n) ")
            if "n" in yn or "N" in yn:
                return False
            if "y" in yn or "Y" in yn:
                rmtree(name)
                mkdir(name)
                return True
    else:
        mkdir(name)
        return True
    
def verify_images(fnlist):
    """Interactively plots all of the images in 'fnlist' to verify their
    contents are correct. Allows user to mark bad images and change image
    header keywords to correct issues.
    
    Inputs:
    fnlist -> List of filenames to be plotted
    
    Outputs:
    newfnlist -> New list of filenames with bad files removed from the list
    
    """

    #Change pyplot shortcut keys
    oldsavekey = plt.rcParams["keymap.save"]
    plt.rcParams["keymap.save"] = ""
    
    #How many total images are there?
    num_images = len(fnlist)
    flagged = np.zeros(len(fnlist),dtype="bool")
    
    #Interactively look at the images
    current_num = 0
    while True:
        image = openfits(fnlist[current_num],mode="update")
        image[0].header["fpverify"]="True"
        verify = Verify_Images_Plot(image[1].data,
                             image[0].header["OBJECT"],
                             flagged[current_num],
                             current_num+1,
                             num_images)
        image.close()
        if verify.key == "a": current_num = current_num-1
        if verify.key == "d": current_num = current_num+1
        if verify.key == "w": flagged[current_num] = True
        if verify.key == "s": flagged[current_num] = False
        
        #Breakout condition
        if current_num == -1:
            while True:
                yn = raw_input("Finished flagging images? (y/n) ")
                if "n" in yn or "N" in yn:
                    current_num = 0
                    break
                elif "y" in yn or "Y" in yn:
                    break
        elif current_num == num_images:
            while True:
                yn = raw_input("Finished flagging images? (y/n) ")
                if "n" in yn or "N" in yn:
                    current_num = current_num-1
                    break
                elif "y" in yn or "Y" in yn:
                    break
        if current_num == -1 or current_num == num_images or verify.key == "q": break
    
    #Review flagged images
    print "Reviewing "+str(np.sum(flagged))+" flagged images..."
    newfnlist = []
    for i in range(len(fnlist)):
        if flagged[i]:
            image = openfits(fnlist[i],mode="update")
            review = Review_Images_Plot(image[1].data,image[0].header["OBJECT"],fnlist[i])
            if review.key == "h":
                #Correct header
                while True:
                    newhead = raw_input("New object type for header (ARC, FLAT, etc.): ")
                    yn = raw_input("New object header: '"+newhead+"'. Is this okay? (y/n) ")
                    if ("y" in yn or "Y" in yn) and not ("n" in yn or "N" in yn):
                        image[0].header["OBJECT"] = newhead
                        break
                newfnlist.append(fnlist[i])
            if review.key == "a":
                newfnlist.append(fnlist[i])
            if review.key == "d":
                print fnlist[i]+" removed from file list."
                image[0].header["FPGOOD"] = "False"
            image.close()
        else: newfnlist.append(fnlist[i])
    
    plt.rcParams["keymap.save"] = oldsavekey

    return newfnlist
    
def separate_lists(fnlist):
    """Separates a list of images into separate lists based on their header
    keywords 'OBJECT' and 'FILTER'
    
    Inputs:
    fnlist -> A list containing strings, each the location of a fits image
    
    Outputs:
    flatlist -> List of flatfield images
    list_of_objects -> List of distinct objects
    objectlists -> List of lists, each one all of the files for a particular object
    list_of_filters -> List of distinct filters
    filterlists -> List of lists, each one all the files for a given filter
    
    """
    
    #Open each file and look for flats and new filters/objects
    list_of_objects = []
    list_of_filters = []
    flatlist = []
    for i in range(len(fnlist)):
        image = openfits(fnlist[i])
        #Break if image was previously flagged as bad
        fpgood = image[0].header.get("FPGOOD")
        if not (fpgood is None):
            if "False" in fpgood:
                continue
        #Case for flats
        if image[0].header["OBJECT"]=="FLAT": flatlist.append(fnlist[i])
        #Case for ARCs
        elif image[0].header["OBJECT"]=="ARC":
            if not (image[0].header["FILTER"] in list_of_filters): list_of_filters.append(image[0].header["FILTER"])
        #Case for Objects
        else:
            if not (image[0].header["FILTER"] in list_of_filters): list_of_filters.append(image[0].header["FILTER"])
            if not (image[0].header["OBJECT"] in list_of_objects): list_of_objects.append(image[0].header["OBJECT"])
        image.close()
    list_of_filters = np.array(list_of_filters)
    list_of_objects = np.array(list_of_objects)
        
    #Reopen the files and populate the object lists
    objectlists = []
    filterlists = []
    for i in range(len(list_of_objects)): objectlists.append([])
    for i in range(len(list_of_filters)): filterlists.append([])
    for i in range(len(fnlist)):
        image = openfits(fnlist[i])
        #Break if image was previously flagged as bad
        fpgood = image[0].header.get("FPGOOD")
        if not (fpgood is None):
            if "False" in fpgood:
                continue
        if image[0].header["OBJECT"]!="FLAT":
            if image[0].header["OBJECT"]!="ARC": objectlists[np.where(image[0].header["OBJECT"]==list_of_objects)[0][0]].append(fnlist[i])
            filterlists[np.where(image[0].header["FILTER"]==list_of_filters)[0][0]].append(fnlist[i])
        image.close()
    return flatlist, list_of_objects, objectlists, list_of_filters, filterlists
    
def get_aperture(path):
    """Opens a fits image and has the user click on points near the image
    aperture to fit a circle to those points. Once at least three points
    on the aperture boundary have been clicked, begins displaying a best-fit
    circle to those points overlayed on the image.
    
    Inputs:
    path -> Path to the fits file to open
    
    Outputs:
    axcen -> The aperture x center
    aycen -> The aperture y center
    arad -> The aperture radius
    
    """

    image = openfits(path)
    xs, ys = [], []
    axcen, aycen, arad = image[1].data.shape[1]/2, image[1].data.shape[0]/2, 0
    apguess = np.array([image[1].data.shape[1]/2,image[1].data.shape[0]/2,min(image[1].data.shape[0],image[1].data.shape[1])/2],dtype='float64')
    while True:
        click = Click_Near_Aperture_Plot(image[1].data,xs,ys,axcen,aycen,arad)
        if not (click.xcoo is None) and click.xcoo>26 and click.ycoo>26 and click.xcoo<image[1].data.shape[1]-26 and click.ycoo<image[1].data.shape[0]-26:
            newclick = Zoom_Aperture_Plot(image[1].data[click.ycoo-25:click.ycoo+25,click.xcoo-25:click.xcoo+25])
            if not (newclick.xcoo is None):
                xs.append(newclick.xcoo+click.xcoo-25)
                ys.append(newclick.ycoo+click.ycoo-25)
                if len(xs)>2:
                    apparams = minimize(lambda x: np.sum(((np.array(xs)-x[0])**2+(np.array(ys)-x[1])**2-x[2]**2)**2), apguess).x
                    axcen, aycen, arad = apparams[0], apparams[1], apparams[2]
        elif len(xs)>2 and click.xcoo is None: break
    
    image.close()
    
    return axcen, aycen, arad
    
# def cosmic_ray_remove(fnlist):
#     """Runs the LACosmicx routine on a series of images in order to correct
#     cosmic rays. Creates the header keyword 'fpcosmic' and sets its value to
#     'True'
#     
#     Inputs:
#     fnlist -> A list containing the paths to fits images.
#     
#     """
#     
#     for i in range(len(fnlist)):
#         print "Cosmic-ray correcting image "+str(i+1)+" of "+str(len(fnlist))+": "+fnlist[i]
#         image = openfits(fnlist[i],mode="update")
#         image[0].header["fpcosmic"] = "True"
#         _mask, image[0].data = lacosmicx(image[0].data,verbose=False,cleantype='idw')
#         image[0].data[np.isnan(image[0].data)]=0
#         image.close()
#     
#     return
    
def measure_fwhm(fnlist):
    """Gets the coordinates for a handful of stars from the user interactively,
    then attempts to measure the FWHM of those stars to obtain a value for the
    typical PSF of those images. The FWHMs are appended to the image headers
    as the keyword 'FPFWHM'. This routine assumes the PSFs are circularly
    symmetric Gaussians (which is not always the case for SALT, especially if
    the mirror alignment fell apart during a track!).
    
    Inputs:
    fnlist -> A list of strings, each the location of a fits image. The images
              should all have roughly the same coordinate system, because the
              same stellar coordinates are used for all images in the list.
    
    """

    #Open the images
    imagelist = []
    for i in range(len(fnlist)):
        imagelist.append(openfits(fnlist[i],mode="update"))
    
    xcen = imagelist[0][0].header["fpaxcen"]
    ycen = imagelist[0][0].header["fpaycen"]
    arad = imagelist[0][0].header["fparad"]
    data = imagelist[0][1].data[ycen-arad:ycen+arad,xcen-arad:xcen+arad]
    xgrid, ygrid = np.meshgrid( np.arange(50), np.arange(50) )
    starx, stary, starfwhm = [], [], []
    while True:
        plot = Mark_Stars_Plot(data,starx,stary)
        if plot.key == "q" and len(starx) != 0: break
        elif not (plot.xcoo is None) and plot.key == "d" and len(starx) != 0:
            #Find the star nearest to the coordinates
            distarray = (np.array(starx)-plot.xcoo)**2+(np.array(stary)-plot.ycoo)**2
            index = np.where(distarray == np.min(distarray))[0][0]
            starx.pop(index)
            stary.pop(index)
            starfwhm.pop(index)
        elif not (plot.xcoo is None) and plot.key == "z":
            zoom = Zoom_Mark_Stars_Plot(data[plot.ycoo-25:plot.ycoo+25,plot.xcoo-25:plot.xcoo+25])
            if not (zoom.xcoo is None):
                #Measure the FWHM of that star in every image
                fwhms = np.zeros(len(fnlist))
                for i in range(len(fnlist)):
                    zoomdata = imagelist[i][1].data[plot.ycoo+zoom.ycoo+ycen-arad-50:plot.ycoo+zoom.ycoo+ycen-arad,plot.xcoo+zoom.xcoo+xcen-arad-50:plot.xcoo+zoom.xcoo+xcen-arad]
                    xpeak = xgrid[zoomdata==np.max(zoomdata)][0]
                    ypeak = ygrid[zoomdata==np.max(zoomdata)][0]
                    xdata = zoomdata[ypeak,:]
                    ydata = zoomdata[:,xpeak]
                    xfit = GaussFit(np.arange(len(xdata)),xdata)
                    yfit = GaussFit(np.arange(len(ydata)),ydata)
                    if xfit[3]<50 and yfit[3]<50 and xfit[3]>0 and yfit[3]>0:
                        fwhms[i]=(xfit[3]/2+yfit[3]/2)*np.sqrt(8*np.log(2))
                        
                #Plot the FWHM across the track and ask user if it's okay
                if len(stary)>1:
                    oldfwhms = np.median(np.array(starfwhm), axis=0)
                    oldfwhmserr = np.std(np.array(starfwhm), axis=0)
                elif len(stary)==1:
                    oldfwhms = np.median(np.array(starfwhm), axis=0)
                    oldfwhmserr = None
                else:
                    oldfwhms = None
                    oldfwhmserr = None
                fwhmsplot = FWHMs_Plot(fwhms,oldfwhms,oldfwhmserr)
                if not(fwhmsplot.button is None):
                    #If it's okay, append the coordinates and FWHMs to the lists
                    starx.append(zoom.xcoo+plot.xcoo-25)
                    stary.append(zoom.ycoo+plot.ycoo-25)
                    starfwhm.append(fwhms)
    
    #Append the measured FWHMs to the image headers
    imagefwhm = np.median(np.array(starfwhm), axis=0)
    for i in range(len(fnlist)):
        imagelist[i][0].header["fpfwhm"] = imagefwhm[i]
    
    #Close the images
    for i in range(len(fnlist)):
        imagelist[i].close()
    
    return 
    
def make_median(fnlist,outfile):
    """Creates a median image from a series of images.
    
    Inputs:
    fnlist -> List of strings pointing to fits images
    outfile -> Location of resultant median image
    
    """
    
    imagelist = []
    datalist = []
    for i in range(len(fnlist)):
        imagelist.append(openfits(fnlist[i]))
        datalist.append(imagelist[i][1].data)
    datalist = np.array(datalist)
    
    meddata = np.median(datalist,axis=0)
    
    medhdu = PrimaryHDU(meddata)
    medhdu.writeto(outfile,clobber=True)
    
    for i in range(len(fnlist)):
        imagelist[i].close()
    
    return
    
def pipeline(rawdir = "raw", mode = "halpha"):
    """Runs successive steps of the saltfp data reduction, checking along the
    way to see if each step was successful. This is the main driver program of
    the SALT Fabry-Perot pipeline.
    
    Inputs:
    rawdir -> String, containing the path to the 'raw' directory. By
                  default, this is 'raw'
    mode -> Mode for velocity fitting. Currently the only option is H-Alpha
                line fitting.
                  
    """
    
    #Set rest wave based on the mode called
    if mode == "halpha": rest_wave = 6562.81
    
    #Create product directory
    if isdir("product"):
        while True:
            yn = raw_input("Product directory already exists. Recreate it? (y/n) ")
            if "n" in yn or "N" in yn: break
            elif "y" in yn or "Y" in yn:
                rmtree("product")
                break
    if not isdir("product"):
        #Acquire the list of filenames from the raw directory
        fnlist = sorted(listdir(rawdir))
        for i in range(len(fnlist)):
            fnlist[i] = join(rawdir,fnlist[i])
        #Run the first two steps of imred on the first image
        iraf.pysalt(_doprint=0)
        iraf.saltred(_doprint=0)
        iraf.saltprepare(fnlist[0], "temp.fits", "", createvar=False, badpixelimage="", clobber=True, logfile="temp.log", verbose=True)
        iraf.saltbias("temp.fits", "temp.fits", "", subover=True, trim=True, subbias=False, masterbias="", median=False, function="polynomial", order=5,rej_lo=3.0,rej_hi=5.0, niter=10, plotover=False, turbo=False, clobber=True, logfile="temp.log", verbose=True)
        #Create the bad pixel mask
        image = openfits("temp.fits")
        for i in range(1,len(image)):
            mask = image[i].data!=image[i].data
            image[i].data = 1*mask
        image.writeto("badpixmask.fits",clobber="True")
        image.close()
        #Remove temporary files
        remove("temp.fits")
        remove("temp.log")
        #Run the raw images through the first few data reduction pipeline steps
        imred(fnlist, "product", bpmfile = "badpixmask.fits")
        #Delete the temporary bad pixel mask
        remove("badpixmask.fits")
        #Move these raw images into the product directory
        mkdir("product")
        fnlist = sorted(listdir("."))
        for i in range(len(fnlist)):
            if "mfxgbpP" in fnlist[i]:
                move(fnlist[i], join("product",fnlist[i]))
    #List of files in the product directory
    fnlist = sorted(listdir("product"))
    for i in range(len(fnlist)): fnlist[i] = join("product",fnlist[i])
    
    #No longer dealing with single-extension fits files
#     #Convert the images to single-extension fits files
#     # This creates a new directory from scratch in order to preserve the
#     # original raw directory.
#     print "Converting images to single extension fits images..."
#     productlist = []
#     singextlist = []
#     for i in range(len(fnlist)):
#         if splitext(fnlist[i])[1] == ".fits":
#             productlist.append(join(rawdir,fnlist[i]))
#             singextlist.append(join("singext","s"+fnlist[i]))
#     toggle = create_directory("singext")
#     if toggle:
#         make_single_extension(productlist,singextlist)
#     else: print "Skipping creation of single-extension fits images."
    
    #Manual verification of fits images and headers 
    firstimage = openfits(fnlist[0])
    verify = firstimage[0].header.get("fpverify")
    firstimage.close()
    if verify is None:
        while True:
            yn = raw_input("Manually verify image contents? (Recommended) (y/n) ")
            if "n" in yn or "N" in yn:
                print "Skipping manual verification of image contents (Not recommended)"
                break
            if "y" in yn or "Y" in yn:
                fnlist = verify_images(fnlist)
                break
        
    #Make separate lists of the different fits files
    flatlist, list_of_objs, objlists, list_of_filts, filtlists = separate_lists(fnlist)

    #Masking of pixels outside the aperture
    firstimage = openfits(objlists[0][0])
    axcen = firstimage[0].header.get("fpaxcen")
    firstimage.close()
    if axcen is None:
        print "Masking pixels outside the RSS aperture..."
        axcen, aycen, arad = get_aperture(objlists[0][0])
        aperture_mask(fnlist,axcen,aycen,arad)
    else: print "Images have already been aperture-masked."
    
    #CR Removal is now done in an earlier step (but slower... talk to Steve about this)
#     #Cosmic ray removal
#     firstimage = openfits(objlists[0][0])
#     crtoggle = firstimage[0].header.get("fpcosmic")
#     firstimage.close()
#     if crtoggle is None:
#         print "Cosmic-ray correcting images..."
#         cosmic_ray_remove(singextlist)
#     else: print "Images have already been cosmic-ray corrected."
    
    #Uncertainty images are now an extension
#     #Create uncertainty images for the object files
#     uncertlists = []
#     for i in range(len(objlists)):
#         uncertlists.append([])
#         for j in range(len(objlists[i])):
#             uncertlists[i].append(splitext(objlists[i][j])[0]+'_unc'+splitext(objlists[i][j])[1])
#     if not isfile(uncertlists[0][0]):
#         print "Generating uncertainty images..."
#         for i in range(len(objlists)):
#             for j in range(len(objlists[i])):
#                 print "Writing uncertainty image "+uncertlists[i][j]
#                 image = openfits(objlists[i][j])
#                 image[0].data = np.sqrt(image[0].data)
#                 image.writeto(uncertlists[i][j])
#                 image.close()
#     else: print "Uncertainty images already exist."
    
    #Image Flattening
    firstimage = openfits(objlists[0][0])
    flattog = firstimage[0].header.get("fpflat")
    firstimage.close()
    if flattog is None:
        print "Flattening images..."
        if len(flatlist) == 0:
            while True:
                print "Uh oh! No flatfield exposure found!"
                flatpath = raw_input("Enter path to external flat image: (leave blank to skip flattening) ")
                if flatpath == "" or isfile(flatpath): break
        else:
            combine_flat(flatlist,"flat.fits")
            flatpath = "flat.fits"
        if flatpath != "":
            notflatlist = []
            for objlist in objlists: notflatlist+=objlist
            flatten(notflatlist,flatpath)
        else: print "Skipping image flattening. (Not recommended!)"
    else: print "Images have already been flattened."
    
    #Measure stellar FWHMs
    firstimage = openfits(objlists[0][0])
    fwhm = firstimage[0].header.get("fpfwhm")
    firstimage.close()
    if fwhm is None: dofwhm = True
    else:
        while True:
            yn = raw_input("Seeing FWHM has already been measured. Redo this? (y/n) ")
            if "n" in yn or "N" in yn:
                dofwhm = False
                break
            elif "y" in yn or "Y" in yn:
                dofwhm = True
                break
    if dofwhm:
        print "Measuring seeing FWHMs..."
        for objlist in objlists:
            measure_fwhm(objlist)
    
    #Find image centers using ghost pairs
    for i in range(len(objlists)):
        firstimage = openfits(objlists[i][0])
        xcen = firstimage[0].header.get("fpxcen")
        deghosted = firstimage[0].header.get("fpghost")
        firstimage.close()
        if deghosted is None:
            if xcen is None:
                ghosttog = True
            else:
                while True:
                    yn = raw_input("Optical centers already measured for object "+list_of_objs[i]+". Redo this? (y/n) ")
                    if "n" in yn or "N" in yn:
                        ghosttog = False
                        break
                    elif "y" in yn or "Y" in yn:
                        ghosttog = True
                        break
            if ghosttog:
                print "Identifying optical centers for object "+list_of_objs[i]+". This may take a while for crowded fields..."
                find_ghost_centers(objlists[i])
    
    #Deghost images
    firstimage = openfits(objlists[0][0])
    deghosted = firstimage[0].header.get("fpghost")
    firstimage.close()
    if deghosted is None:
        print "Deghosting images..."
        for i in range(len(objlists)):
            for j in range(len(objlists[i])):
                deghost(objlists[i][j],g=0.04)
    else: print "Images have already been deghosted."

    #Make separate directories for each object.
    # This is the first bit since 'singext' to create a new directory, because
    # this is the first point where it's really necessary to start treating the
    # images from different tracks very differently.
    for i in range(len(objlists)):
        if isdir(list_of_objs[i].replace(" ", "")):
            while True:
                yn = raw_input("A directory for object "+list_of_objs[i]+" already exists. Recreate? (y/n) ")
                if "n" in yn or "N" in yn:
                    do_copy = False
                    break
                elif "y" in yn or "Y" in yn:
                    do_copy = True
                    rmtree(list_of_objs[i].replace(" ", ""))
                    break
        else: do_copy = True
        if do_copy:
            mkdir(list_of_objs[i].replace(" ", ""))
            for j in range(len(objlists[i])):
                copyfile(objlists[i][j],join(list_of_objs[i].replace(" ", ""),split(objlists[i][j])[1]))
        for j in range(len(objlists[i])):
            objlists[i][j] = join(list_of_objs[i].replace(" ", ""),split(objlists[i][j])[1])
    #Update the filter lists (THIS IS TERRIBLE AND I SHOULD HAVE USED POINTERS AND I'M SO SORRY)
    for i in range(len(filtlists)):
        for j in range(len(filtlists[i])):
            for k in range(len(objlists)):
                for l in range(len(objlists[k])):
                    if split(filtlists[i][j])[1]==split(objlists[k][l])[1]: filtlists[i][j] = objlists[k][l]
    
    #Image alignment and normalization
    for i in range(len(objlists)):
        firstimage = openfits(objlists[i][0])
        aligned = firstimage[0].header.get("fpphot")
        firstimage.close()
        if aligned is None:
            print "Aligning and normalizing images for object "+list_of_objs[i]+"..."
            align_norm(objlists[i])
        else: print "Images for object "+list_of_objs[i]+" have already been aligned and normalized."
        
    #Make a median image for each object
    for i in range(len(objlists)):
        if isfile(join(list_of_objs[i].replace(" ", ""),"median.fits")):
            while True:
                yn = raw_input("Median image for object "+list_of_objs[i]+" already exists. Replace it? (y/n) ")
                if "n" in yn or "N" in yn:
                    break
                elif "y" in yn or "Y" in yn:
                    make_median(objlists[i],join(list_of_objs[i].replace(" ", ""),"median.fits"))
                    break
        else: make_median(objlists[i],join(list_of_objs[i].replace(" ", ""),"median.fits"))

    #Wavelength calibrations
    for i in range(len(list_of_filts)):
        #Do a separate calibration for each filter
        firstimage = openfits(filtlists[i][0])
        calf = firstimage[0].header.get("fpcalf")
        firstimage.close()
        if not(calf is None):
            while True:
                yn = raw_input("Wavelength solution already found for filter "+list_of_filts[i]+". Redo it? (y/n) ")
                if "n" in yn or "N" in yn:
                    break
                elif "y" in yn or "Y" in yn:
                    fit_wave_soln(filtlists[i])
                    break
        else: fit_wave_soln(filtlists[i])
    
    #Sky ring removal
    for i in range(len(objlists)):
        #Check to see if sky rings have already been removed
        firstimage = openfits(objlists[i][0])
        deringed = firstimage[0].header.get("fpdering")
        firstimage.close()
        if deringed is None:
            print "Subtracting sky rings for object "+list_of_objs[i]+"..."
            sub_sky_rings(objlists[i],[join(list_of_objs[i].replace(" ", ""),"median.fits")]*len(objlists[i]))
        else: print "Sky ring subtraction already done for object "+list_of_objs[i]
    
    #Creation of data cube and convolution to uniform PSF
    for i in range(len(objlists)):
        if isdir(list_of_objs[i].replace(" ", "")+"_cube"):
            while True:
                yn = raw_input("A data cube for object "+list_of_objs[i]+" already exists. Recreate? (y/n) ")
                if "n" in yn or "N" in yn:
                    do_create = False
                    break
                elif "y" in yn or "Y" in yn:
                    do_create = True
                    rmtree(list_of_objs[i].replace(" ", "")+"_cube")
                    break
        else: do_create = True
        if do_create:
            mkdir(list_of_objs[i].replace(" ", "")+"_cube")
            for j in range(len(objlists[i])):
                image = openfits(objlists[i][j])
                fwhm = image[0].header["fpfwhm"]
                if j==0: largestfwhm = fwhm
                if fwhm>largestfwhm: largestfwhm=fwhm
                image.close()
            while True:
                user_fwhm = raw_input("Enter desired final fwhm or leave blank to use default ("+str(largestfwhm)+" pix) ")
                if user_fwhm == "":
                    user_fwhm = largestfwhm
                    break
                else:
                    try: user_fwhm = float(user_fwhm)
                    except ValueError: print "That wasn't a valid number..."
                    else:
                        if user_fwhm<largestfwhm:
                            print "Final fwhm must exceed "+str(largestfwhm)+" pixels."
                        else: break
            desired_fwhm = user_fwhm
            for j in range(len(objlists[i])):
                make_final_image(objlists[i][j], #Input image
                                 join(list_of_objs[i].replace(" ", "")+"_cube",split(objlists[i][j])[1]), #Output image
                                 desired_fwhm, #Desired final FWHM
                                 clobber=True)
    
    #Get final lists for the velocity map fitting for each object
    final_lists = []
    for i in range(len(list_of_objs)):
        final_lists.append([])
        for j in range(len(objlists[i])):
            final_lists[i].append(join(list_of_objs[i].replace(" ", "")+"_cube",split(objlists[i][j])[1]))
            
    #Shift to solar velocity frame
    for i in range(len(list_of_objs)):
        firstimage = openfits(final_lists[i][0])
        velshift = firstimage[0].header.get("fpsolar")
        firstimage.close()
        if velshift is None:
            print "Performing solar velocity shift for object "+list_of_objs[i]+"..."
            solar_velocity_shift(final_lists[i], rest_wave)
        else: print "Solar velocity shift for object "+list_of_objs[i]+" already done."
    
    #Velocity map fitting
    for i in range(len(list_of_objs)):
        if (isfile(join(list_of_objs[i].replace(" ", "")+"_cube","velocity.fits"))):
            while True:
                yn = raw_input("Velocity map already fitted for object "+list_of_objs[i]+". Redo this? (y/n) ")
                if "n" in yn or "N" in yn:
                    domap = False
                    break
                elif "y" in yn or "Y" in yn:
                    domap = True
                    break
        else: domap = True
    if domap:
        print "Fitting velocity map for object "+list_of_objs[i]+"..."
        if mode == "halpha":
            fit_velmap_ha_n2_mode(final_lists[i],
                                  list_of_objs[i].replace(" ", "")+"_cube",
                                  clobber=True)
    
    #Clean velocity map
    for i in range(len(list_of_objs)):
        make_clean_map(list_of_objs[i].replace(" ", "")+"_cube", clobber=True)

if __name__ == "__main__":
    if len(sys.argv)==2: pipeline(rawdir=sys.argv[1], mode="halpha")
    else: pipeline(mode="halpha")
