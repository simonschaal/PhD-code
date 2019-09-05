import numpy as np
from lmfit import Model


# Fit functions
# LORENTZ: offset, amplitude, resonance, FWHM
def lorentz(x, A0, A, res, w):
    return A0 + 2*A/np.pi * w/(w**2+4*(x-res)**2)

# FANO: offset, amplitude, FWHM, resonance
def fano(x, A0, A, q, gamma, res):
    return A0 + A*(q*gamma/2+x-res)**2/((gamma/2)**2+(x-res)**2)

# Convert logarithmic VNA data to linear refection coefficient
def SxxtoGamma(y, bgmin=0, bgmax=-1):
    ylin=10**(y/20) # divided by 20 because reflection coefficient is defined in terms of voltages
    ylin=ylin/np.mean(ylin[bgmin:bgmax]) # subtract BG based on selection
    return ylin

# perform a fit using lmfit given inital values and bounds
def fit(x, y, func, init={}, bounds={}):
    
    gmodel = Model(func)
    params = gmodel.make_params()
    
    for param in init.keys():
        if param in params:
            params[param].value = init[param]

    for param in bounds.keys():
        if param in params:
            params[param].min = init[param][0]
            params[param].max = init[param][1]
            
    return gmodel.fit(y, params, x=x)

# calcualte Q from resonance fit results based on FWHM and resonance frequency
def getQ(fitresults, keys=dict(freq='res', fwhm='w')):
    return fitresults.best_values[keys['freq']]/fitresults.best_values[keys['fwhm']]
    