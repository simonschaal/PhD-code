import numpy as np

# calculate probabislity desnity for singlet and triplet
def probST(V, signal, noise, tau, pt=0.5 , t1=4.5e-3):#4.5
    nS=(1-pt)/np.sqrt(2*np.pi*noise**2) * np.exp( -(V+signal/2)**2/(2*noise**2) )
    intrg=np.arange(-signal/2, signal/2, signal/1000)
    integrand=lambda x: tau/t1 * pt/signal * 1/np.sqrt(2*np.pi*noise**2) * np.exp( -(x[:,np.newaxis]+signal/2)/(signal) * tau/t1) * np.exp( -(V-x[:,np.newaxis])**2/(2*noise**2) )
    integral=np.trapz(integrand(intrg), intrg, axis=0)
    nT=np.exp(-tau/t1)/np.sqrt(2*np.pi*noise**2) * pt * np.exp( -(V-signal/2)**2/(2*noise**2) ) + integral
    
    return nS, nT

# calculate fidelity from probability densities
def fidelity(Vt, signal, noise, tau, t1=4.5e-3):
    # fs
    integrand=lambda x: probST(x, signal, noise, tau, t1=t1)[0]
    intrg=np.arange(Vt, 3*signal, signal/1000) #10signal
    fS=1-np.trapz(integrand(intrg), intrg)
    
    # fT
    integrand=lambda x: probST(x, signal, noise, tau, t1=t1)[1]
    intrg=np.arange(-3*signal, Vt, signal/1000)
    fT=1-np.trapz(integrand(intrg), intrg)
    
    return fS, fT

# find optimal threshold
def optimizeVt(Vtmin, Vtmax, steps, signal, noise, tau, t1):
    Vts=np.linspace(Vtmin, Vtmax, steps)
    fmean=np.zeros(len(Vts))
    for i in range(len(Vts)):
        fmean[i]=np.mean(fidelity(Vts[i], signal, noise, tau, t1=t1))
        
    imax=np.argmax(fmean)
    return fmean, Vts[imax], fmean[imax]