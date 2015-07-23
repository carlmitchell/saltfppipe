import numpy as np
import matplotlib.pyplot as plt
from astropy.io.fits import open as openfits
from os.path import join, split, isfile
from sys import exit as crash
from gauss_fit import GaussFit, gaussfunc
from saltfppipe.wavelib import get_libraries
from fpfunc import fpfunc_for_curve_fit
from scipy.optimize.minpack import curve_fit
from saltfppipe.fpfunc import fpfunc_for_curve_fit_with_t

class WaveSolnPlot:
    def __init__(self,r,z,t,wave,res,colors,time_divs):
        
        rms = np.sqrt(np.average(res**2))
        
        self.key = None
        self.xcoo = None
        self.ycoo = None
        self.axis = None
        
        fig = plt.figure()
        
        ax1 = fig.add_subplot(221)
        self.plot1 = ax1.scatter(z,res,color=colors)
        ax1.set_xlabel('z')
        ax1.set_ylabel('Residual')
        
        ax2 = fig.add_subplot(222)
        self.plot2 = ax2.scatter(r,res,color=colors)
        ax2.set_xlabel('Radius [pix]')
        ax2.set_ylabel('Residual')
        
        ax3 = fig.add_subplot(223)
        self.plot3 = ax3.scatter(t,res,color=colors)
        ax3.set_xlabel('Time [min]')
        ax3.set_ylabel('Residual')
        for i in range(len(time_divs)):
            ax3.axvline(time_divs[i],color='black')
        
        ax4 = fig.add_subplot(224)
        self.plot4 = ax4.scatter(wave,res,color=colors)
        ax4.set_xlabel('Wavelength [Ang]')
        ax4.set_ylabel('Residual')
        
        plt.suptitle("Red points are ARC rings. Blue points are sky rings.\n"+
                     "Press 'd' to delete nearest point, 'r' to restore all points.\n"+
                     "Press 't' to toggle linear time dependence.\n"+
                     "Press 'w' to add a piecewise time break, 'q' to remove all time breaks.\n"+
                     "Press 'a' to accept the solution and move on.\n"+
                     "RMS = "+str(rms))
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.75)
        
        self.cid1 = self.plot1.figure.canvas.mpl_connect('key_press_event',self.keypress)
        self.cid2 = self.plot2.figure.canvas.mpl_connect('key_press_event',self.keypress)
        self.cid3 = self.plot3.figure.canvas.mpl_connect('key_press_event',self.keypress)
        self.cid4 = self.plot4.figure.canvas.mpl_connect('key_press_event',self.keypress)
    
        plt.show()
        
    def keypress(self,event):
        if event.key in ["d","t","q","a","r","w"]:
            self.key = event.key
            if event.inaxes in [self.plot1.axes,self.plot2.axes,self.plot3.axes,self.plot4.axes]:
                self.xcoo = event.xdata
                self.ycoo = event.ydata
                if event.inaxes == self.plot1.axes: self.axis = 1
                if event.inaxes == self.plot2.axes: self.axis = 2
                if event.inaxes == self.plot3.axes: self.axis = 3
                if event.inaxes == self.plot4.axes: self.axis = 4
            plt.close()
            return
        
class PlotRingProfile:
    def __init__(self,data,rbins,intbins,xcen,ycen,radii,numstring):
        self.key = None
        self.xcoo = None
        self.ycoo = None
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        self.plot = ax1.imshow(data,
                               cmap = "Greys",
                               aspect = "equal",
                               vmin = np.percentile(data,5),
                               vmax = np.percentile(data,95),
                               origin = "lower",
                               interpolation = "nearest")
        ax1.plot(rbins,intbins,color="blue")
        for i in range(len(radii)):
            ax1.plot(xcen+radii[i]*np.cos(np.linspace(0,2*np.pi,100)),ycen+radii[i]*np.sin(np.linspace(0,2*np.pi,100)),color="red")
        ax1.set_xlim([0,data.shape[1]])
        ax1.set_ylim([0,data.shape[0]])
        plt.title("Image "+numstring+"\n"+
                  "Press 'w' to fit a ring ('e' to force-mark a ring)\n"+
                  "Press 's' to remove the nearest marked ring.\n"+
                  "Press 'a'/'d' to change images.")
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
    
    def keypress(self,event):
        if event.key in ["w","e","a","s","d"]:
            self.key = event.key
            if event.inaxes == self.plot.axes:
                self.xcoo = event.xdata
                self.ycoo = event.ydata
            plt.close()
            return

class PlotRingFit:
    def __init__(self,x,y,fit):
        self.key = None
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        self.plot = ax1.scatter(x,y,color="blue")
        ax1.plot(np.linspace(np.min(x),np.max(x),1000),gaussfunc(np.linspace(np.min(x),np.max(x),1000), fit[0], fit[1], fit[2], fit[3]), color="red")
        plt.axvline(fit[2],color='red')
        plt.title("Press 'w' to accept ring fit.\n"+
                  "Press 's' to reject it.")
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
        
    def keypress(self,event):
        if event.key in ["w","s"]:
            self.key = event.key
            plt.close()
            return
        
class TimePlot:
    def __init__(self,t,res,colors,time_divs):
        self.xcoo = None
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        self.plot = ax1.scatter(t,res,color=colors)
        for i in range(len(time_divs)):
            ax1.axvline(time_divs[i],color="black")
        plt.title("Black lines show current time divisions.\n"+
                  "Press 'w' where you would like to add a time division.\n"
                  "Press 'q' to not add a time division.")
        ax1.set_xlabel("Time [min]")
        ax1.set_ylabel("Residual [Ang]")
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
        
    def keypress(self,event):
        if event.key in ["w","q"]:
            if event.key == "w" and event.inaxes == self.plot.axes:
                self.xcoo = event.xdata
            plt.close()
            return

def fit_wave_soln(fnlist):
    """Fits a wavelength solution to rings in a set of images. Appends this
    wavelength solution to the image headers as the keywords:
    'fpcala', 'fpcalb', ... 'fpcalf'
    
    Each object in fnlist must have a corresponding "median.fits" in its
    image directory, or this routine will not work.
    
    ARC ring images are fitted by adjusting the center, while the center is held
    fixed for night sky rings. A combination of fits to both sets of rings is
    used to determine a wavelength solution for the whole set of images.
    
    If the ARC rings disagree substantially with the night sky rings, it is
    recommended that users delete the ARC rings from the fit and use only the
    night sky rings.
    
    It is also known that the wavelength solution can sometimes be piecewise in
    time when a large jump in 'z' happens between two images; i.e. different
    wavelength solutions exist before and after the jump. The routine allows
    the user to make a piecewise solution for this reason, but this ability
    should be used sparingly.
    
    This routine contains one of the few hard-coded numbers in the pipeline,
    Fguess=5600. Currently F values are not written to the fits image headers,
    and this is a reasonable guess.
    
    Inputs:
    fnlist -> List of strings, each the path to a fits image. These images
    should all have been taken with the same order filter. If not, the routine
    will crash.
    
    """
    
    #This bit takes care of the 's' to save shortcut in matplotlib.
    oldsavekey = plt.rcParams["keymap.save"]
    plt.rcParams["keymap.save"] = ""
        
    #Open all of the images
    imagelist = []
    arclist = []
    objlist = []
    for i in range(len(fnlist)):
        imagelist.append(openfits(fnlist[i]))
        if i == 0: filt = imagelist[0][0].header["FILTER"]
        if imagelist[i][0].header["FILTER"] != filt:
            print "Error! Some of these images are in different filters!"
            crash()
        if imagelist[i][0].header["OBJECT"]=="ARC": arclist.append(imagelist[i])
        else:
            if not isfile(join(split(fnlist[i])[0],"median.fits")):
                print "Error! No 'median.fits' file found."
                crash()
            medimage = openfits(join(split(fnlist[i])[0],"median.fits"))
            imagelist[i][0].data += -medimage[0].data
            medimage.close()
            objlist.append(imagelist[i])
    
    #Load wavelength libraries
    arclib, nightlib = get_libraries(filt)
    if arclib is None:
        print "Error! Your filter isn't the wavelength library!"
        crash()

    



    #This next bit fits all of the rings that the user marks

    #Fit rings in the object images
    radlists = []
    for i in range(len(objlist)):
        radlists.append([])
    i=0
    while True:
        xgrid, ygrid = np.meshgrid(np.arange(objlist[i][0].data.shape[1]), np.arange(objlist[i][0].data.shape[0]))
        xcen = objlist[i][0].header["FPXCEN"]
        ycen = objlist[i][0].header["FPYCEN"]
        axcen = objlist[i][0].header["FPAXCEN"]
        aycen = objlist[i][0].header["FPAYCEN"]
        arad = objlist[i][0].header["FPARAD"]
        rgrid = np.sqrt((xgrid - xcen)**2 + (ygrid - ycen)**2)
        rbins = np.arange(arad-np.int(max(abs(axcen-xcen),abs(aycen-ycen))))+1
        intbins = np.empty_like(rbins)
        for j in range(len(rbins)):
            intbins[j] = np.median(objlist[i][0].data[np.logical_and(np.logical_and(objlist[i][0].data!=0,rgrid<rbins[j]),rgrid>rbins[j]-1)])
        ringplot = PlotRingProfile(objlist[i][0].data[aycen-arad:aycen+arad,axcen-arad:axcen+arad], #Data to be plotted. Only want stuff inside aperture
                                   rbins+(xcen-axcen)+arad, #Radii bins shifted to image center
                                   intbins*arad/np.percentile(np.abs(intbins),98)+(ycen-aycen)+arad, #Intensity bins, rescaled and shifted by image center
                                   xcen-axcen+arad, ycen-aycen+arad, #Shifted center
                                   radlists[i], #Previously fitted rings
                                   repr(i+1)+"/"+repr(len(objlist))) #numstring
        #Changing images and loop breakout conditions
        if ringplot.key == "d": i+=1
        if ringplot.key == "a": i+=-1
        if i == -1 or i == len(objlist):
            while True:
                yn = raw_input("Finished marking sky rings? (y/n) ")
                if "n" in yn or "N" in yn:
                    if i == -1: i=0
                    if i == len(objlist): i = len(objlist)-1
                    break
                elif "y" in yn or "Y" in yn:
                    break
        if i == -1 or i == len(objlist): break
        #Force-marking a ring
        if ringplot.key == "e" and ringplot.xcoo != None: radlists[i].append(ringplot.xcoo-arad-(xcen-axcen))
        #Deleting a ring
        if ringplot.key == "s" and ringplot.xcoo != None and len(radlists[i])>0:
            radlists[i].pop(np.argmin(np.array(radlists[i])-np.sqrt((ringplot.xcoo-arad-(xcen-axcen))**2 + (ringplot.ycoo-arad-(ycen-aycen))**2)))
        #Fitting a ring profile
        if ringplot.key == "w" and ringplot.xcoo != None:
            x = rbins[max(ringplot.xcoo-arad-(xcen-axcen)-50,0):min(ringplot.xcoo-arad-(xcen-axcen)+50,len(rbins))]**2
            y = intbins[max(ringplot.xcoo-arad-(xcen-axcen)-50,0):min(ringplot.xcoo-arad-(xcen-axcen)+50,len(rbins))]
            fit = GaussFit(x,y)
            fitplot = PlotRingFit(x,y,fit)
            if fitplot.key == "w": radlists[i].append(np.sqrt(fit[2]))
    zo = []
    to = []
    ro = []
    for i in range(len(objlist)):
        for j in range(len(radlists[i])):
            zo.append(objlist[i][0].header["ET1Z"])
            to.append(objlist[i][0].header["JD"])
            ro.append(radlists[i][j])
            
    #Fit rings in the ARC images
    xcen = objlist[0][0].header["FPXCEN"]
    ycen = objlist[0][0].header["FPYCEN"]
    radlists = []
    for i in range(len(arclist)):
        radlists.append([])
    i=0
    while True:
        xgrid, ygrid = np.meshgrid(np.arange(arclist[i][0].data.shape[1]), np.arange(arclist[i][0].data.shape[0]))
        axcen = arclist[i][0].header["FPAXCEN"]
        aycen = arclist[i][0].header["FPAYCEN"]
        arad = arclist[i][0].header["FPARAD"]
        rgrid = np.sqrt((xgrid - xcen)**2 + (ygrid - ycen)**2)
        rbins = np.arange(arad-np.int(max(abs(axcen-xcen),abs(aycen-ycen))))+1
        intbins = np.empty_like(rbins)
        for j in range(len(rbins)):
            intbins[j] = np.median(arclist[i][0].data[np.logical_and(np.logical_and(arclist[i][0].data!=0,rgrid<rbins[j]),rgrid>rbins[j]-1)])
        ringplot = PlotRingProfile(arclist[i][0].data[aycen-arad:aycen+arad,axcen-arad:axcen+arad], #Data to be plotted. Only want stuff inside aperture
                                   rbins+(xcen-axcen)+arad, #Radii bins shifted to image center
                                   intbins*arad/np.percentile(np.abs(intbins),98)+(ycen-aycen)+arad, #Intensity bins, rescaled and shifted by image center
                                   xcen-axcen+arad, ycen-aycen+arad, #Shifted center
                                   radlists[i], #Previously fitted rings
                                   repr(i+1)+"/"+repr(len(arclist))) #numstring
        #Changing images and loop breakout conditions
        if ringplot.key == "d": i+=1
        if ringplot.key == "a": i+=-1
        if i == -1 or i == len(arclist):
            while True:
                yn = raw_input("Finished marking ARC rings? (y/n) ")
                if "n" in yn or "N" in yn:
                    if i == -1: i=0
                    if i == len(arclist): i = len(arclist)-1
                    break
                elif "y" in yn or "Y" in yn:
                    break
        if i == -1 or i == len(arclist): break
        #Force-marking a ring
        if ringplot.key == "e" and ringplot.xcoo != None: radlists[i].append(ringplot.xcoo-arad-(xcen-axcen))
        #Deleting a ring
        if ringplot.key == "s" and ringplot.xcoo != None and len(radlists[i])>0:
            radlists[i].pop(np.argmin(np.array(radlists[i])-np.sqrt((ringplot.xcoo-arad-(xcen-axcen))**2 + (ringplot.ycoo-arad-(ycen-aycen))**2)))
        #Fitting a ring profile
        if ringplot.key == "w" and ringplot.xcoo != None:
            x = rbins[max(ringplot.xcoo-arad-(xcen-axcen)-50,0):min(ringplot.xcoo-arad-(xcen-axcen)+50,len(rbins))]**2
            y = intbins[max(ringplot.xcoo-arad-(xcen-axcen)-50,0):min(ringplot.xcoo-arad-(xcen-axcen)+50,len(rbins))]
            fit = GaussFit(x,y)
            fitplot = PlotRingFit(x,y,fit)
            if fitplot.key == "w": radlists[i].append(np.sqrt(fit[2]))
    za = []
    ta = []
    ra = []
    for i in range(len(arclist)):
        for j in range(len(radlists[i])):
            za.append(arclist[i][0].header["ET1Z"])
            ta.append(arclist[i][0].header["JD"])
            ra.append(radlists[i][j])
    
#     #Load previous ring fits from a text file - COMMENT THIS OUT LATER
#     rr,zz,tt = np.loadtxt("test.out",unpack=True)
#     za = list(zz[zz>0])
#     ta = list(tt[zz>0])
#     ra = list(rr[zz>0])
#     zo = list(zz[zz<0])
#     to = list(tt[zz<0])
#     ro = list(rr[zz<0])

    #Now we try to get a good guess at the wavelengths
    
    #Get a good guess at which wavelengths are which
    Bguess = objlist[0][0].header["ET1B"]
    Fguess = 5600
    
    #Figure out A by matching rings to the wavelength libraries
    master_r = np.array(ro+ra)
    master_z = np.array(zo+za)
    wavematch = np.zeros_like(master_r)
    isnight = np.array([True]*len(ro)+[False]*len(ra))
    oldrms = 10000 #Really high initial RMS for comparisons
    for i in range(len(master_r)):
        if isnight[i]: lib = nightlib
        else: lib = arclib
        for j in range(len(lib)):
            #Assume the i'th ring is the j'th line
            Aguess = lib[j]*np.sqrt(1+master_r[i]**2/Fguess**2)-Bguess*master_z[i]
            #What are all of the other rings, given this A?
            waveguess = (Aguess+Bguess*master_z)/np.sqrt(1+master_r**2/Fguess**2)
            for k in range(len(master_r)):
                if isnight[k]: wavematch[k] = nightlib[np.argmin(np.abs(nightlib-waveguess[k]))]
                else: wavematch[k] = arclib[np.argmin(np.abs(arclib-waveguess[k]))]
            rms = np.sqrt(np.average((waveguess-wavematch)**2))
            if rms < oldrms:
                #This is the new best solution. Keep it!
                oldrms = rms
                bestA = Aguess
                master_wave = wavematch.copy()
    
    #Make more master arrays for the plotting
    master_t = np.array(to+ta)
    t0 = np.min(master_t)
    master_t += -t0
    master_t *= 24*60 #Convert to minutes
    master_color = np.array(len(ro)*["blue"]+len(ra)*["red"]) #Colors for plotting
    toggle = np.ones(len(master_r),dtype="bool")
    dotime = False
    time_dividers = []
    
    #Do the interactive plotting
    while True:
        rplot = master_r[toggle]
        zplot = master_z[toggle]
        tplot = master_t[toggle]
        colorplot = master_color[toggle]
        waveplot = master_wave[toggle]
        fitplot = np.zeros(len(waveplot))
        xs = np.zeros((3,len(rplot)))
        xs[0] = rplot
        xs[1] = zplot
        xs[2] = tplot
        fit = [0]*(len(time_dividers)+1)
        time_dividers = sorted(time_dividers)
        if len(time_dividers)>1: print "Warning: Too many time divisions is likely unphysical. Be careful!"
        for i in range(len(time_dividers)+1):
            #Create a slice for all of the wavelengths before this time divider
            #but after the one before it
            if len(time_dividers)==0: tslice = tplot==tplot
            elif i == 0: tslice = tplot<time_dividers[i]
            elif i==len(time_dividers): tslice = tplot>time_dividers[i-1]
            else: tslice = np.logical_and(tplot<time_dividers[i],tplot>time_dividers[i-1])
            if dotime:
                fit[i] = curve_fit(fpfunc_for_curve_fit_with_t, xs[:,tslice], waveplot[tslice], p0=(bestA,Bguess,0,Fguess))[0]
                fitplot[tslice] = fpfunc_for_curve_fit_with_t(xs[:,tslice], fit[i][0], fit[i][1], fit[i][2], fit[i][3])
            else:
                fit[i] = curve_fit(fpfunc_for_curve_fit, xs[:,tslice], waveplot[tslice], p0=(bestA,Bguess,Fguess))[0]
                fitplot[tslice] = fpfunc_for_curve_fit(xs[:,tslice], fit[i][0], fit[i][1], fit[i][2])
        resid = waveplot - fitplot
        solnplot = WaveSolnPlot(rplot,zplot,tplot,waveplot,resid,colorplot,time_dividers)
        #Breakout case
        if solnplot.key == "a":
            while True:
                for i in range(len(time_dividers)+1):
                    if dotime: print "Solution 1: A = "+str(fit[i][0])+", B = "+str(fit[i][1])+", E = "+str(fit[i][2])+", F = "+str(fit[i][3])
                    else: print "Solution 1: A = "+str(fit[i][0])+", B = "+str(fit[i][1])+", F = "+str(fit[i][2])
                print "Residual rms="+str(np.sqrt(np.average(resid**2)))+" for "+repr(len(time_dividers)+1)+" independent "+repr(3+dotime)+"-parameter fits to "+repr(len(rplot))+" rings."
                yn = raw_input("Accept wavelength solution? (y/n) ")
                if "n" in yn or "N" in yn:
                    break
                elif "y" in yn or "Y" in yn:
                    solnplot.key = "QUIT"
                    break
        if solnplot.key == "QUIT": break
        #Restore all points case
        if solnplot.key == "r": toggle = np.ones(len(master_r),dtype="bool")
        #Delete nearest point case
        if solnplot.key == "d" and solnplot.axis != None:
            #Figure out which plot was clicked in
            if solnplot.axis == 1:
                #Resid vs. z plot
                z_loc = solnplot.xcoo
                resid_loc = solnplot.ycoo
                dist2 = ((zplot-z_loc)/(np.max(zplot)-np.min(zplot)))**2 + ((resid-resid_loc)/(np.max(resid)-np.min(resid)))**2
            elif solnplot.axis == 2:
                #Resid vs. R plot
                r_loc = solnplot.xcoo
                resid_loc = solnplot.ycoo
                dist2 = ((rplot-r_loc)/(np.max(rplot)-np.min(rplot)))**2 + ((resid-resid_loc)/(np.max(resid)-np.min(resid)))**2
            elif solnplot.axis == 3:
                #Resit vs. T plot
                t_loc = solnplot.xcoo
                resid_loc = solnplot.ycoo
                dist2 = ((tplot-t_loc)/(np.max(tplot)-np.min(tplot)))**2 + ((resid-resid_loc)/(np.max(resid)-np.min(resid)))**2
            elif solnplot.axis == 4:
                #Resid vs. Wave plot
                wave_loc = solnplot.xcoo
                resid_loc = solnplot.ycoo
                dist2 = ((waveplot-wave_loc)/(np.max(waveplot)-np.min(waveplot)))**2 + ((resid-resid_loc)/(np.max(resid)-np.min(resid)))**2
            #Get the radius and time of the worst ring
            r_mask = rplot[dist2 == np.min(dist2)][0]
            t_mask = tplot[dist2 == np.min(dist2)][0]
            toggle[np.logical_and(master_r == r_mask, master_t == t_mask)] = False
        #Fit for time case
        if solnplot.key == "t": dotime = not dotime
        #Add time break
        if solnplot.key == "w":
            timeplot = TimePlot(tplot,resid,colorplot,time_dividers)
            if timeplot.xcoo != None: time_dividers.append(timeplot.xcoo)
        #Remove time breaks
        if solnplot.key == "q":
            time_dividers = []

    #Close all images
    for i in range(len(fnlist)):
        imagelist[i].close()
        
    #For each image, write the central wavelength and F to the image header
    for i in range(len(fnlist)):
        image = openfits(fnlist[i],mode="update")
        image_t = (image[0].header["JD"]-t0)*24*60
        #Figure out which time division it's in
        div_index = np.where(np.array(time_dividers)>image_t)[0]
        if len(div_index>0): div_index = div_index[0]
        else: div_index = len(time_dividers)
        image_fit = fit[div_index]
        if dotime:
            image_wave0 = image_fit[0]+image_fit[1]*image[0].header["ET1Z"]+image_fit[2]*image_t
            image_F = image_fit[3]
        else:
            image_wave0 = image_fit[0]+image_fit[1]*image[0].header["ET1Z"]
            image_F = image_fit[2]
        image[0].header["FPWAVE0"] = image_wave0
        image[0].header["FPCALF"] = image_F
        image.close()
    
    #Restore the old keyword shortcut
    plt.rcParams["keymap.save"] = oldsavekey
    
    
    
    return