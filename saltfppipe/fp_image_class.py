from astropy.io import fits
import numpy as np
from sys import exit as crash

class FPImage:
    """Class to wrap various operations for a SALT FP image file, in particular
    provides a nice wrapper for the 3 or 4 image extensions and some useful
    header keywords.
    
    Inputs:
    filename -> String. Path to the fits file to be opened.
    update -> Bool. Whether to overwrite changes to original file. Default=False.
    
    """
    
    def __init__(self, filename, update=False):
        
        #Open the file
        self.openimage = fits.open(filename)
        
        #Store its own filename and update parameters
        self.filename = filename
        self.update = update
        
        #Assign data arrays
        self.inty = self.openimage[1].data
        self.vari = self.openimage[2].data
        self.badp = self.openimage[3].data
        try: self.wave = self.openimage[4].data
        except: self.wave = None
        
        #Xsize and Ysize of image
        self.xsize = self.inty.shape[1]
        self.ysize = self.inty.shape[0]
        
        #Assign header
        self.header = self.openimage[0].header
        
        #Header keywords from the telescope
        if "Etalon 1" in self.header.get("et-state"): self.state = "mr"
        elif "Etalon 2" in self.header.get("et-state"): self.state = "lr"
        #else: crash("Etalon mode [kw: ET-STATE] could not be determined for"+
        #            self.filename)
        else: self.state = 'mr'
        self.filter = self.header.get("filter")
        self.object = self.header.get("object")
        self.jd = self.header.get("jd")
        if self.state == "mr":
            self.z = self.header.get("et1z")
            self.a = self.header.get("et1a")
            self.b = self.header.get("et1b")
            self.f = self.header.get("et1f")
        elif self.state == "lr":
            self.z = self.header.get("et2z")
            self.a = self.header.get("et2a")
            self.b = self.header.get("et2b")
            self.f = self.header.get("et2f")
        self.bins = self.header.get("ccdsum")
        self.datestring = self.header.get("date-obs")
        self.ut = self.header.get("utc-obs")
        self.ra = self.header.get("ra")
        self.dec = self.header.get("dec")
        
        #Non-boolean header keywords from the pipeline
        self.xcen = self.header.get("fpxcen")
        self.ycen = self.header.get("fpycen")
        self.axcen = self.header.get("fpaxcen")
        self.aycen = self.header.get("fpaycen")
        self.arad = self.header.get("fparad")
        self.wave0 = self.header.get("fpwave0")
        self.fwhm = self.header.get("fpfwhm")
        self.calf = self.header.get("fpcalf")
        self.solarvel = self.header.get("fpsolar")
        self.dnorm = self.header.get("fpdnorm")
        self.calrms = self.header.get("fpcalrms")
        self.calnring = self.header.get("ncalring")
        self.calnpars = self.header.get("ncalpars")
        self.calnfits = self.header.get("ncalfits")

        self.xshift = self.header.get("fpxshift")
        self.yshift = self.header.get("fpyshift")
        
        #Boolean header keywords
        self.phottog = self.header.get("fpphot")
        self.ghosttog = self.header.get("fpghost")
        self.flattog = self.header.get("fpflat")
        self.ringtog = self.header.get("fpdering")
        self.verifytog = self.header.get("fpverify")
        self.goodtog = self.header.get("fpgood")
        
        #Rescale F by the pixel binning
        if not self.f is None: self.f /= int(self.bins.split()[0])
        
        #If F is zero, it shouldn't be. Someone just effed up the headers.
        if self.f == 0: self.f = None
        
    def close(self):
        """Closes the open file. If the file was opened in update mode,
        overwrites the old file.

        """
        
        if self.update==True: self.writeto(self.filename, clobber=True)
        self.openimage.close()
        
    def writeto(self, outfilepath, clobber=False):
        """Writes the open FP image file to the specified output path. If
        the output path already exists, overwrites the file only if clobber
        is True. (default False). Updates non-None header values.
        
        """
        
        def update_header_kw(header, kwstring, newval):
            """Method for updating header values if they are not None."""
            if not newval is None: header[kwstring] = newval
            return header
        
        #Update telescope headers
        update_header_kw(self.header, "filter", self.filter)
        update_header_kw(self.header, "object", self.object)
        update_header_kw(self.header, "jd", self.jd)
        update_header_kw(self.header, "ccdsum", self.bins)
        update_header_kw(self.header, "date-obs", self.datestring)
        update_header_kw(self.header, "utc-obs", self.ut)
        update_header_kw(self.header, "ra", self.ra)
        update_header_kw(self.header, "dec", self.dec)
        if self.state == "mr":
            update_header_kw(self.header, "et1z", self.z)
            update_header_kw(self.header, "et1a", self.a)
            update_header_kw(self.header, "et1b", self.b)
            if self.f is None:
                update_header_kw(self.header, "et1f", self.f)
            else:
                update_header_kw(self.header, "et1f",
                                 self.f*int(self.bins.split()[0]))
        elif self.state == "lr":
            update_header_kw(self.header, "et2z", self.z)
            update_header_kw(self.header, "et2a", self.a)
            update_header_kw(self.header, "et2b", self.b)
            if self.f is None:
                update_header_kw(self.header, "et2f", self.f)
            else:
                update_header_kw(self.header, "et2f",
                                 self.f*int(self.bins.split()[0]))
        
        #Update non-boolean pipeline headers
        update_header_kw(self.header, "fpxcen", self.xcen)
        update_header_kw(self.header, "fpycen", self.ycen)
        update_header_kw(self.header, "fpaxcen", self.axcen)
        update_header_kw(self.header, "fpaycen", self.aycen)
        update_header_kw(self.header, "fparad", self.arad)
        update_header_kw(self.header, "fpwave0", self.wave0)
        update_header_kw(self.header, "fpfwhm", self.fwhm)
        update_header_kw(self.header, "fpcalf", self.calf)
        update_header_kw(self.header, "fpsolar", self.solarvel)
        update_header_kw(self.header, "fpdnorm", self.dnorm)
        update_header_kw(self.header, "fpcalrms", self.calrms)
        update_header_kw(self.header, "ncalpars", self.calnpars)
        update_header_kw(self.header, "ncalring", self.calnring)
        update_header_kw(self.header, "ncalfits", self.calnfits)
        update_header_kw(self.header, "fpxshift", self.xshift)
        update_header_kw(self.header, "fpyshift", self.yshift)
        
        #Update boolean pipeline headers
        update_header_kw(self.header, "fpphot", self.phottog)
        update_header_kw(self.header, "fpghost", self.ghosttog)
        update_header_kw(self.header, "fpflat", self.flattog)
        update_header_kw(self.header, "fpdering", self.ringtog)
        update_header_kw(self.header, "fpverify", self.verifytog)
        update_header_kw(self.header, "fpgood", self.goodtog)
        
        #Update the actual image data itself
        self.openimage[0].header = self.header
        self.openimage[1].data = self.inty
        self.openimage[2].data = self.vari
        self.openimage[3].data = self.badp
        if not self.wave is None: self.openimage[4].data = self.wave
        
        #Write to the output file
        self.openimage.writeto(outfilepath,clobber=clobber)

    def xyarrays(self):
        """Returns arrays of X and Y coordinates in the image. I.e. for a 2x3
        image will return:
        
        [[0, 1, 2], and [[0, 0, 0],
         [0, 1, 2]]      [1, 1, 1]]
        
        """
        
        return np.meshgrid(np.arange(self.xsize),
                           np.arange(self.ysize))
        
    def rarray(self, xcen, ycen):
        """Returns array of radius values at about a specified point. I.e. for
        a 3x3 array with specified center (xcen,ycen) = (1,1), will return:
        
        [[1.412, 1.000, 1.412],
         [1.000, 0.000, 1.000],
         [1.412, 1.000, 1.412]]
        
        """
        
        return np.sqrt((self.xyarrays()[0]-xcen)**2+
                       (self.xyarrays()[1]-ycen)**2)
        
    def skybackground(self):
        """Fits for the sky background level and standard deviation in the
        intensity array of the image, then returns these values. Does so by
        computing the mean and st. dev. of all of the pixels in the 10th to
        90th percentile of the image.
        
        Outputs:
        skyavg, skysig, skyvar
        
        """
        
        data = self.inty[self.badp==0]
        lowval = np.percentile(data,10)
        highval = np.percentile(data,90)
        data = data[np.logical_and(data>lowval,data<highval)]
        
        return np.average(data), np.std(data), np.std(data)/np.sqrt(len(data))
    
        