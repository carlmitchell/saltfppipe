import numpy as np

def get_libraries(filt):
    if filt == "PI06530": return neon_6530(), night_6530()
    elif filt == "PI06645": return neon_6645(), night_6645()
    
    else: return None, None
    
def neon_6530():
    return np.array([6402.248,
                     6506.5281,
                     6532.8822,
                     6598.9529
                     ])
    
def night_6530():
    return np.array([6464.824, #There are about 5 unresolved OH- lines here, the two brightest are 6465.847 and 6463.801. This is their average.
                     6470.974,
                     6477.910,
                     6498.729,
                     6504.995,
                     6522.425,
                     6533.044,
                     6544.022,
                     6553.617,
                     6562.82, #H-alpha
                     6568.779,
                     6577.284, #Average of an OH- doublet
                     6583.41, #[NII]
                     6596.643,
                     6604.134 #Average of an OH- doublet
                     ])
    
def neon_6645():
    return np.array([6532.8822,
                     6598.9529,
                     6652.0927,
                     6678.2762,
                     6717.0430
                     ])
    
def night_6645():
    return np.array([6544.022,
                     6553.617,
                     6562.82, #H-alpha
                     6568.779,
                     6577.284, #Average of an OH- doublet
                     6583.41, #[NII]
                     6596.643,
                     6604.134, #Average of an OH- doublet
                     6627.624, #Very dim line, unlikely to find it
                     6634.230, #Average of an OH- doublet
                     6661.771, #Average of an OH- doublet
                     6667.627, #Average of an OH- doublet
                     6699.156, #Average of an OH- doublet
                     6704.385 #Average of an OH- doublet
                     ])