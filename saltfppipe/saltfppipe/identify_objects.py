from pyraf import iraf
from os import remove

def identify_objects(fnlist,skysiglist,fwhmlist,suffix=".coo"):
    """
    DEPRECIATED. THIS ROUTINE NOW DOES NOTHING.
    
    Runs the IRAF routine 'daofind' to locate objects in a series of images,
    creating coordinate files.
    
    Inputs:
    fnlist -> List of strings, each the path to a fits image.
    skysiglist -> List of floats, each the sky background sigma for an image.
    fwhmlist -> List of floats, each the FWHM of objects in an image.
    suffix -> Suffix for the coordinate files. '.coo' by default.
    
    Outputs:
    coolist -> List of strings, each the path to the coordinate files created.
    
    """
    
    print "Identifying objects in images..."
    
    coolist = []
    
    #Open IRAF packages
    iraf.noao(_doprint=0)
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    for i in range(len(fnlist)):
        coolist.append(fnlist[i]+suffix)
        iraf.daofind(image=fnlist[i],
                     output=fnlist[i]+suffix,
                     fwhmpsf=fwhmlist[i],
                     sigma=skysiglist[i],
                     threshold=4.0,
                     datamin='INDEF',
                     datamax='INDEF',
                     verify='N')
    return coolist
    
def clean_files(fnlist,suffix=".coo"):
    """
    DEPRECIATED. THIS ROUTINE NOW DOES NOTHING.

    Deletes the coordinate files produced by identify_objects.
    
    Inputs:
    fnlist -> List of strings, each the path to a fits image.
    suffix -> Suffix for the coordinate files. '.coo' by default.
    
    """
    
    print "Deleting temporary coordinate files..."
    
    for i in range(len(fnlist)):
        remove(fnlist[i]+suffix)
        
    return