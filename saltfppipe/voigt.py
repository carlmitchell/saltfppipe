import numpy as np
from scipy.optimize.minpack import curve_fit


def nVoigt(x, wave0, sigma, gamma):
    """Normalized Voigt Profile at the point 'x', centered around 'wave0',
    with gaussian width 'sigma' and lorentzian width 'gamma'.

    Magic numbers come from:
    P. Thompson, D.E. Cox, J.B. Hastings, J. Appl. Cryst. 1987, 20, 79.

    """

    # Calculate the new full-width half-max (magic numbers come from paper in
    # documentation)
    f = fwhm(sigma, gamma)

    # Calculate the mixing parameter (magic numbers come from paper in
    # documentation)
    m = mix(gamma, f)

    # Calculate normalized Gaussian and Lorentzians with FWHM 'f' at the point
    # of interest
    g = (2./f)*np.sqrt(np.log(2.)/np.pi)*np.exp(-4.*np.power((x-wave0)/f, 2))
    l = (f/(2.*np.pi))/(np.power(x-wave0, 2)+np.power(f/2., 2))

    # Create the normalized Voigt profile at the point of interest
    v = m*l + (1.-m)*g

    return v


def fwhm(sigma_g, gamma_l):
    """Calculates the new FWHM of a pseudo-Voigt profile according to the work
    of P. Thompson, D.E. Cox, J.B. Hastings, J. Appl. Cryst. 1987, 20, 79.

    Inputs:
    sigma_g -> Gaussian standard deviation sigma
    gamma_l -> Lorentzian FWHM gamma

    Outputs:
    f -> The new fwhm

    """

    s = 2.*np.sqrt(2.*np.log(2.))
    return np.power(np.power(s*sigma_g, 5)*np.power(gamma_l, 0) +
                    np.power(s*sigma_g, 4)*np.power(gamma_l, 1)*2.69296 +
                    np.power(s*sigma_g, 3)*np.power(gamma_l, 2)*2.42843 +
                    np.power(s*sigma_g, 2)*np.power(gamma_l, 3)*4.47163 +
                    np.power(s*sigma_g, 1)*np.power(gamma_l, 4)*0.07842 +
                    np.power(s*sigma_g, 0)*np.power(gamma_l, 5), 0.2)


def mix(gamma_l, f):
    """Calculates the mixing parameter 'm' for mixing Gaussian and Lorentzian
    components together.

    Magic numbers come from:
    P. Thompson, D.E. Cox, J.B. Hastings, J. Appl. Cryst. 1987, 20, 79.

    Inputs:
    gamma_l -> The lorentzian FWHM gamma
    f -> The desired total FWHM

    Outputs:
    m-> The mixing parameter

    """

    return (1.36603*np.power(gamma_l/f, 1) -
            0.47719*np.power(gamma_l/f, 2) +
            0.11116*np.power(gamma_l/f, 3))


def d2Voigt(prof, x):
    """Evaluates the second derivatives of a Voigt profile V(C,I,l_0,s_g,g_l)
    at the location x.

    The profile is actually a pseudo-voigt, based on the work of
    P. Thompson, D.E. Cox, J.B. Hastings, J. Appl. Cryst. 1987, 20, 79.

    Inputs:
    prof -> [Cont, Int, lambda_0, sigma_g, gamma_l]
    x -> point at which to evaluate the second derivatives

    Outputs:
    d2V -> 5x5 array containing the Hessian matrix [[d2V/dC2, d2V/dCdI, etc.]]

    """

    [_C, I, lambda_0, sigma_g, gamma_l] = prof

    # Calculate the new FWHM and its first and second derivatives:
    f = fwhm(sigma_g, gamma_l)
    df_dsigma = (0.2*2.*np.sqrt(2.*np.log(2.))/np.power(f,4)) * ( 5.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,4)*np.power(gamma_l,0) + 
                                                                  4.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,3)*np.power(gamma_l,1)*2.69296 + 
                                                                  3.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,2)*np.power(gamma_l,2)*2.42843 + 
                                                                  2.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,1)*np.power(gamma_l,3)*4.47163 + 
                                                                  1.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,0)*np.power(gamma_l,4)*0.07842 )
    df_dgamma = (0.2/np.power(f,4)) * ( 1.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,4)*np.power(gamma_l,0)*2.69296 + 
                                        2.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,3)*np.power(gamma_l,1)*2.42843 + 
                                        3.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,2)*np.power(gamma_l,2)*4.47163 + 
                                        4.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,1)*np.power(gamma_l,3)*0.07842 + 
                                        5.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,0)*np.power(gamma_l,4) )
    d2f_dsigma2 = ( (0.2*np.power(2.*np.sqrt(2.*np.log(2.)),2)/np.power(f,4)) * ( 5.*4.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,3)*np.power(gamma_l,0) + 
                                                                                  4.*3.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,2)*np.power(gamma_l,1)*2.69296 + 
                                                                                  3.*2.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,1)*np.power(gamma_l,2)*2.42843 + 
                                                                                  2.*1.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,0)*np.power(gamma_l,3)*4.47163 ) - 
                                                                                4.*np.power(df_dsigma,2)/f )
    d2f_dgamma2 = ( (0.2/np.power(f,4)) * ( 2.*1.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,3)*np.power(gamma_l,0)*2.42843 + 
                                            3.*2.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,2)*np.power(gamma_l,1)*4.47163 + 
                                            4.*3.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,1)*np.power(gamma_l,2)*0.07842 + 
                                            5.*4.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,0)*np.power(gamma_l,3) ) - 
                                          4.*np.power(df_dgamma,2)/f )                            
    d2f_dsigma_dgamma = ( (0.2*2.*np.sqrt(2.*np.log(2.))/np.power(f,4)) * ( 4.*1.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,3)*np.power(gamma_l,0)*2.69296 + 
                                                                          3.*2.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,2)*np.power(gamma_l,1)*2.42843 + 
                                                                          2.*3.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,1)*np.power(gamma_l,2)*4.47163 + 
                                                                          1.*4.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,0)*np.power(gamma_l,3)*0.07842 ) - 
                                                                          4.*df_dgamma*df_dsigma/f )

    # Calculate the mixing parameter and its first and second PARTIAL
    # derivatives:
    m = mix(gamma_l, f)
    dm_dgamma = (1.*1.36603*np.power(gamma_l, 0)/np.power(f, 1) -
                 2.*0.47719*np.power(gamma_l, 1)/np.power(f, 2) +
                 3.*0.11116*np.power(gamma_l, 2)/np.power(f, 3))
    dm_df = (-1.*1.36603*np.power(gamma_l, 1)/np.power(f, 2) +
             2.*0.47719*np.power(gamma_l, 2)/np.power(f, 3) -
             3.*0.11116*np.power(gamma_l, 3)/np.power(f, 4))
    d2m_df2 = (1.*2.*1.36603*np.power(gamma_l, 1)/np.power(f, 3) -
               2.*3.*0.47719*np.power(gamma_l, 2)/np.power(f, 4) +
               3.*4.*0.11116*np.power(gamma_l, 3)/np.power(f, 5))
    d2m_dgamma_df = (-1.*1.*1.36603*np.power(gamma_l, 0)/np.power(f, 2) +
                     2.*2.*0.47719*np.power(gamma_l, 1)/np.power(f, 3) -
                     3.*3.*0.11116*np.power(gamma_l, 2)/np.power(f, 4))

    # Calculate the normalized Gaussian and Lorentzian profiles at the
    # location of interest, as well as their derivatives
    g = ((4./f)*np.sqrt(np.log(2.)/np.pi) *
         np.exp(-4.*np.power((x-lambda_0)/f, 2)))
    l = (f/(2.*np.pi))/(np.power(x-lambda_0, 2)+np.power(f/2., 2))
    dl_dlambda = 4.*np.pi*(x-lambda_0)*np.power(l, 2)/f
    dg_dlambda = 8.*(x-lambda_0)*g/np.power(f, 2)
    dl_df = l/f - np.pi*np.power(l, 2)
    dg_df = (g/f)*(8*np.power((x-lambda_0)/f, 2)-1)

    # Calculate a few derivatives w.r.t. the parameter q = gamma/f
    # (because I apparently can't do math)
    dm_dq = (1*1.36603*np.power(gamma_l/f, 0) -
             2*0.47719*np.power(gamma_l/f, 1) +
             3*0.11116*np.power(gamma_l/f, 2))
    d2m_dq2 = (-1*2*0.47719*np.power(gamma_l/f, 0) +
               2*3*0.11116*np.power(gamma_l/f, 1))
    dq_dgamma = 1/f - gamma_l*df_dgamma/np.power(f, 2)

    # Calculate first derivatives of the Voigt profile at the location of
    # interest
    dV_dlambda = I*(m*dl_dlambda + (1.-m)*dg_dlambda)
    dV_dsigma = I*df_dsigma*(m*dl_df+(l-g)*dm_df+(1.-m)*dg_df)
    dV_dgamma = I*(m*dl_df*df_dgamma+(1.-m)*dg_df*df_dgamma +
                   (l-g)*dm_dq*dq_dgamma)

    # Calculate the second derivatives of the Voigt profile at the location of
    # interest
    d2V_dI_dlambda = dV_dlambda/I
    d2V_dI_dsigma = dV_dsigma/I
    d2V_dI_dgamma = dV_dgamma/I
    d2V_dlambda2 = I*(m*(2*l*dl_dlambda*4*np.pi*(x-lambda_0)/f-4*np.pi*np.power(l,2)/f)+(1.-m)*(8*(x-lambda_0)*dg_dlambda/np.power(f,2)-8*g/np.power(f,2)))
    d2V_dlambda_dsigma = (4*I*(x-lambda_0)/np.power(f,3))*df_dsigma*(dm_df*(np.pi*np.power(l*f,2)-2*g*f)+np.pi*m*(np.power(l,2)*f-2*np.pi*np.power(l,3)*np.power(f,2))+(1-m)*(-4*g+2*g*(8*np.power((x-lambda_0)/f,2)-1)))
    d2V_dlambda_dgamma = (4*I*(x-lambda_0)/np.power(f,3))*(m*np.pi*df_dgamma*(f*np.power(l,2)-2*np.pi*np.power(f,2)*np.power(l,3))+g*df_dgamma*(1.-m)*(16*np.power((x-lambda_0)/f,2)-6)+(dm_df*df_dgamma+dm_dgamma)*(np.pi*np.power(f*l,2)-2*g*f))
    d2V_dsigma2 = I*d2f_dsigma2*(m*dl_df+(1.-m)*dg_df+(l-g)*dm_df)+I*np.power(df_dsigma,2)*(m*(dl_df/f-l/np.power(f,2)-2*np.pi*l*dl_df)+dm_df*dl_df+(l-g)*d2m_df2+dm_df*(dl_df-dg_df)-dm_df*dg_df+(1.-m)*(dg_df/f)*(8*np.power((x-lambda_0)/f,2)-1)-(1.-m)*(g/np.power(f,2))*(8*np.power((x-lambda_0)/f,2)-1)-(1.-m)*(16*g/np.power(f,4))*np.power(x-lambda_0,2))
    d2V_dsigma_dgamma = I*d2f_dsigma_dgamma*(m*dl_df+(1.-m)*dg_df+(l-g)*dm_df)+I*df_dsigma*((dm_dgamma+dm_df*df_dgamma)*(dl_df-dg_df)+(d2m_dgamma_df+d2m_df2*df_dgamma)*(l-g)+df_dgamma*(m*(dl_df/f-l/np.power(f,2)-2*np.pi*l*dl_df)+dm_df*(dl_df-dg_df)+(1.-m)*((dg_df/f)*(8*np.power((x-lambda_0)/f,2)-1)-(g/np.power(f,2))*(8*np.power((x-lambda_0)/f,2)-1)-16*g*np.power(x-lambda_0,2)/np.power(f,4))))
    d2V_dgamma2 = I*(dm_dq*dq_dgamma*dl_df*df_dgamma+m*(dl_df/f-2*np.pi*l*dl_df-l/np.power(f,2))*np.power(df_dgamma,2)+m*dl_df*d2f_dgamma2-dm_dq*dq_dgamma*dg_df*df_dgamma+(1.-m)*(dg_df*(8*np.power((x-lambda_0)/f,2)-1)/f-g*(8*np.power((x-lambda_0)/f,2)-1)/np.power(f,2)-16*g*np.power(x-lambda_0,2)/np.power(f,4))*np.power(df_dgamma,2)+(1.-m)*dg_df*d2f_dgamma2+(dl_df-dg_df)*df_dgamma*dm_dq*dq_dgamma+(l-g)*d2m_dq2*np.power(dq_dgamma,2)+(l-g)*dm_dq*(2*gamma_l*np.power(df_dgamma,2)/np.power(f,3)-2*df_dgamma/np.power(f,2)-gamma_l*d2f_dgamma2/np.power(f,2)))

    d2V = np.asarray([[np.zeros_like(x),np.zeros_like(x),np.zeros_like(x),np.zeros_like(x),np.zeros_like(x)],
                      [np.zeros_like(x),np.zeros_like(x),d2V_dI_dlambda,d2V_dI_dsigma,d2V_dI_dgamma],
                      [np.zeros_like(x),d2V_dI_dlambda,d2V_dlambda2,d2V_dlambda_dsigma,d2V_dlambda_dgamma],
                      [np.zeros_like(x),d2V_dI_dsigma,d2V_dlambda_dsigma,d2V_dsigma2,d2V_dsigma_dgamma],
                      [np.zeros_like(x),d2V_dI_dgamma,d2V_dlambda_dgamma,d2V_dsigma_dgamma,d2V_dgamma2]])

    return d2V


def dVoigt(prof, x, _y=[], _dy=[]):
    """Evaluates the first derivatives of a Voigt profile V(C,I,l_0,s_g,g_l)
    at the location x.

    The profile is actually a pseudo-voigt, based on the work of
    P. Thompson, D.E. Cox, J.B. Hastings, J. Appl. Cryst. 1987, 20, 79.

    Inputs:
    prof -> [Cont, Int, lambda_0, sigma_g, gamma_l]
    x -> point at which to evaluate the derivatives

    _y -> Dummy argument for the minimization routine
    _dy -> Dummy argument for the minimization routine

    Outputs:
    dV -> 5-long array containing [dV/dC, dV/dI, etc.]

    """

    [_C, I, lambda_0, sigma_g, gamma_l] = prof

    # New FWHM and Mixing Parameter
    f = fwhm(sigma_g, gamma_l)
    m = mix(gamma_l, f)

    # Their derivatives
    df_dsigma = (0.2*2.*np.sqrt(2.*np.log(2.))/np.power(f,4)) * ( 5.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,4)*np.power(gamma_l,0) + 
                                                                  4.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,3)*np.power(gamma_l,1)*2.69296 + 
                                                                  3.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,2)*np.power(gamma_l,2)*2.42843 + 
                                                                  2.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,1)*np.power(gamma_l,3)*4.47163 + 
                                                                  1.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,0)*np.power(gamma_l,4)*0.07842 )
    df_dgamma = (0.2/np.power(f,4)) * ( 1.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,4)*np.power(gamma_l,0)*2.69296 + 
                                        2.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,3)*np.power(gamma_l,1)*2.42843 + 
                                        3.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,2)*np.power(gamma_l,2)*4.47163 + 
                                        4.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,1)*np.power(gamma_l,3)*0.07842 + 
                                        5.*np.power(2.*np.sqrt(2.*np.log(2.))*sigma_g,0)*np.power(gamma_l,4) )
    dm_df = (-1.*1.36603*np.power(gamma_l, 1)/np.power(f, 2) +
             2.*0.47719*np.power(gamma_l, 2)/np.power(f, 3) -
             3.*0.11116*np.power(gamma_l, 3)/np.power(f, 4))
    dm_dq = (1*1.36603*np.power(gamma_l/f, 0) -
             2*0.47719*np.power(gamma_l/f, 1) +
             3*0.11116*np.power(gamma_l/f, 2))
    dq_dgamma = 1/f - gamma_l*df_dgamma/np.power(f, 2)

    # Normalized Gaussian and Lorentzian and their derivatives
    g = ((4./f)*np.sqrt(np.log(2.)/np.pi) *
         np.exp(-4.*np.power((x-lambda_0)/f, 2)))
    l = (f/(2.*np.pi))/(np.power(x-lambda_0, 2)+np.power(f/2., 2))
    dl_dlambda = 4.*np.pi*(x-lambda_0)*np.power(l, 2)/f
    dg_dlambda = 8.*(x-lambda_0)*g/np.power(f, 2)
    dl_df = l/f - np.pi*np.power(l, 2)
    dg_df = (g/f)*(8*np.power((x-lambda_0)/f, 2)-1)

    # Normalized Voigt
    v = m*l + (1.-m)*g

    # Various derivatives
    dV_dC = np.ones_like(x)
    dV_dI = v
    dV_dlambda = I*(m*dl_dlambda + (1.-m)*dg_dlambda)
    dV_dsigma = I*df_dsigma*(m*dl_df+(l-g)*dm_df+(1.-m)*dg_df)
    dV_dgamma = (I*(m*dl_df*df_dgamma+(1.-m) *
                    dg_df*df_dgamma+(l-g)*dm_dq*dq_dgamma))

    dV = np.asarray([dV_dC, dV_dI, dV_dlambda, dV_dsigma, dV_dgamma])

    return dV


def Voigt(prof, x):
    """Evaluates a Voigt Profile V(C,I,l_0,s_g,g_l) at the location x.

    The profile is actually a pseudo-voigt, based on the work of
    P. Thompson, D.E. Cox, J.B. Hastings, J. Appl. Cryst. 1987, 20, 79.

    Inputs:
    prof -> [Cont, Int, lambda_0, sigma_g, gamma_l]
    x -> point at which to evaluate the profile

    Outputs:
    V -> Value of the Voigt profile at that point.

    """

    # The various profile parameters
    [C, I, lambda_0, sigma_g, gamma_l] = prof

    # Rescale the normalized profile by the intensity and add on the continuum
    V = C+I*nVoigt(x, lambda_0, sigma_g, gamma_l)

    return V


def chi2(prof, x, y, dy):
    """Calculates the chi-squared statistic for a pseudo-Voigt profile
    Voigt(prof,x) with a set of data y with error bars dy.

    Inputs:
    prof -> Line profile [C, I, l_0, s_g, g_l]
    x -> Wavelength array
    y -> Intensity array
    dy -> Error bars on intensity
    """

    return np.sum(np.power((Voigt(prof, x)-y)/dy, 2))


def FitVoigt(x, y, dy, fix=[None, None, None, None, None]):
    """Fits a single Voigt profile to a set of data points using a maximum-
    likelihood method. Returns the values of the best fitting fit parameters,
    as well as the Fisher uncertainty matrix for those parameters.

    The parameters of the fit are:
    [Continuum, Intensity, Peak Wavelength, Sigma_g, Gamma_l]

    Inputs:
    x -> Array of x coordinates for the data points (wavelengths)
    y -> Array of y coordinates at the x coordinates in 'x' (intensities)
    dy -> Array of uncertainties in the y coordinates in 'y'
    fix -> 5-long array of either 'None' or a numerical value. 'None' means a
           parameter will be fitted for. A value means the parameter will be
           held fixed at that value.
           e.g. fix= [None, None, None, None, 1.5] will fit for the first four
           parameters while holding Gamma_l fixed at a value of 1.5

    Returns:
    fit -> Fitted values for the 5 parameters which were fitted
           (including fixed ones)
    unc -> Uncertainties in the 5 parameters. Uncertainties in fixed parameters
           are zero.
    redchi2 -> Reduced chi^2 statistic chi^2/n, where n is the number of fitted
               parameters

    """
    return


def ha_and_n2_voigts(x, c, i1, i2, waveha, sig1, sig2, gam1, gam2):
    ha = 6562.81
    n2 = 6583.41
    beta = ((waveha/ha)**2-1)/((waveha/ha)**2+1)
    waven2 = n2*(1+beta)/np.sqrt(1-beta**2)
    return c+i1*nVoigt(x, waveha, sig1, gam1)+i2*nVoigt(x, waven2, sig2, gam2)


def fit_ha_and_n2(x, y, dy):
    ha = 6562.81
    n2 = 6583.41
    mask = np.logical_not(np.logical_or(np.isnan(y), np.isnan(x)))
    args = np.argsort(x[mask])
    x = x[mask][args]
    y = y[mask][args]
    if dy is not None:
        dy = dy[mask][args]
    guess = np.zeros(8)
    guess[0] = np.percentile(y, 10)
    guess[3] = x[y == np.max(y)][0]
    guess_beta = ((guess[3]/ha)**2-1)/((guess[3]/ha)**2+1)
    guess_n2 = n2*(1+guess_beta)/np.sqrt(1-guess_beta**2)
    haslice = np.logical_and(x > guess[3]-10, x < guess[3]+10)
    n2slice = np.logical_and(x > guess_n2-10, x < guess_n2+10)
    guess[1] = np.sum(0.5*(y[haslice][1:]+y[haslice][:-1]-2*guess[0]) *
                      np.diff(x[haslice]))
    guess[2] = np.sum(0.5*(y[n2slice][1:]+y[n2slice][:-1]-2*guess[0]) *
                      np.diff(x[n2slice]))
    guess[4:8] = 2
    try:
        fit = curve_fit(ha_and_n2_voigts, x, y, p0=guess, sigma=dy)[0]
        success = True
    except:
        fit = guess
        success = False
    return fit, success


def one_Voigt(x, cont, inty, wave0, sig, gam):
    return cont + inty*nVoigt(x, wave0, sig, gam)


def fit_only_ha(x, y, dy):
    mask = np.logical_not(np.logical_or(np.isnan(y), np.isnan(x)))
    args = np.argsort(x[mask])
    x = x[mask][args]
    y = y[mask][args]
    if dy is not None:
        dy = dy[mask][args]
    guess = np.zeros(5)
    guess[0] = np.percentile(y, 10)
    guess[2] = x[y == np.max(y)][0]
    guess[1] = np.sum(0.5*(y[1:]+y[:-1]-2*guess[0])*np.diff(x))
    guess[3:5] = 2
    try:
        fit = curve_fit(one_Voigt, x, y, p0=guess, sigma=dy)[0]
        success = True
    except:
        fit = guess
        success = False

    return fit, success
