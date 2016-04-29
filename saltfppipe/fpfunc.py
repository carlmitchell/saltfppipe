def fpfunc_for_curve_fit(x, A, B, F):
    """A version of FPFunc for Curve_Fit, where x=[r,z,t]

    """

    return (A+B*x[1])/(1+x[0]**2/F**2)**0.5


def fpfunc_for_curve_fit_with_t(x, A, B, E, F):
    """A version of FPFunc made for Curve_Fit, where x=[r,z,t]

    """

    return (A+B*x[1]+E*x[2])/(1+x[0]**2/F**2)**0.5
