from pyraf import iraf

def make_single_extension(fnlist, newfnlist):
    """Converts fits images to the single-extension format that is more
    compatible with IRAF, CFITSIO, etc. Does this via the pysalt task
    'salt2iraf'.
    
    Inputs:
    fnlist -> List of strings, each one the location of a multi-extension image
    newfnlist -> List of strings, locations for the new single-extension images
    
    """

    #Open various iraf packages
    iraf.pysalt(_doprint=0)
    iraf.saltred(_doprint=0)
    
    #Run salt2iraf on each image
    for i in range(len(fnlist)):
        iraf.salt2iraf(images=fnlist[i],outimages=newfnlist[i],outpref="")
    
    return