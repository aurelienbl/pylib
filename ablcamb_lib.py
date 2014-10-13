

import numpy as np
import matplotlib.pylab as plt

## We use CAMB convention for ordering the Cls, i.e.


## scalCls : TT, EE, TE, PP, PT
## lensedCls : TT, EE, BB, TE
## lenspotentialCls : TT, EE, BB, TE, dd, dT, dE
 

def load_lensedCls(filename, lmax):
    "load_lensedCls(filename, lmax)"

    lensedCls=np.zeros((4, lmax+1))

    dummy=np.loadtxt(filename)

    lensedCls[0,2:lmax+1]=dummy[:,1][0:lmax-1]
    lensedCls[1,2:lmax+1]=dummy[:,2][0:lmax-1]
    lensedCls[2,2:lmax+1]=dummy[:,3][0:lmax-1]
    lensedCls[3,2:lmax+1]=dummy[:,4][0:lmax-1]
   
    l=np.arange(lmax+1)

    ## Rescaling TT
    lensedCls[0,2:lmax+1]= lensedCls[0,2:lmax+1]/(l[2:lmax+1]*(l[2:lmax+1]+1))*(2*np.pi)
    lensedCls[1,2:lmax+1]= lensedCls[1,2:lmax+1]/(l[2:lmax+1]*(l[2:lmax+1]+1))*(2*np.pi)
    lensedCls[2,2:lmax+1]= lensedCls[2,2:lmax+1]/(l[2:lmax+1]*(l[2:lmax+1]+1))*(2*np.pi)
    lensedCls[3,2:lmax+1]= lensedCls[3,2:lmax+1]/(l[2:lmax+1]*(l[2:lmax+1]+1))*(2*np.pi)
   
    return lensedCls



def load_scalCls(filename, lmax):
    "load_scalCls(filename, lmax):"
    scalCls=np.zeros((5, lmax+1))

    dummy=np.loadtxt(filename)

    scalCls[0,2:lmax+1]=dummy[:,1][0:lmax-1]
    scalCls[1,2:lmax+1]=dummy[:,2][0:lmax-1]
    scalCls[2,2:lmax+1]=dummy[:,3][0:lmax-1]
    scalCls[3,2:lmax+1]=dummy[:,4][0:lmax-1]
    scalCls[4,2:lmax+1]=dummy[:,5][0:lmax-1]
    l=np.arange(lmax+1)

    ## Rescaling TT
    scalCls[0,2:lmax+1]= scalCls[0,2:lmax+1]/(l[2:lmax+1]*(l[2:lmax+1]+1))*(2*np.pi)
    scalCls[1,2:lmax+1]= scalCls[1,2:lmax+1]/(l[2:lmax+1]*(l[2:lmax+1]+1))*(2*np.pi)
    scalCls[2,2:lmax+1]= scalCls[2,2:lmax+1]/(l[2:lmax+1]*(l[2:lmax+1]+1))*(2*np.pi)
    scalCls[3,2:lmax+1]= scalCls[3,2:lmax+1]/(l[2:lmax+1]*1.0)**4/ 7.4311e12
    scalCls[4,2:lmax+1]= scalCls[4,2:lmax+1]/(l[2:lmax+1]*1.0)**3/ 7.4311e12

    return scalCls



def save_scalCls(scalCls,filename, lmax):
    "save_scalCls(scalCls,filename, lmax):"
    print 'toto'


    scalCls_for_camb=np.zeros((6, lmax-1))
    lcamb=np.arange(2, lmax+1)
    l=np.arange(lmax+1)
    llpi=(1.0*l)*(l+1)/2.0/np.pi
    scalCls_for_camb[0,:]=lcamb
    scalCls_for_camb[1,:]=(scalCls[0,:]*llpi)[2:]
    scalCls_for_camb[2,:]=(scalCls[1,:]*llpi)[2:]
    scalCls_for_camb[3,:]=(scalCls[2,:]*llpi)[2:]
    scalCls_for_camb[4,:]=(scalCls[3,:]*(l*1.0)**4*7.4311e12)[2:]
    scalCls_for_camb[5,:]=(scalCls[4,:]*(l*1.0)**3*7.4311e12)[2:]

    np.savetxt(filename, scalCls_for_camb.T, fmt=('%d', '%15.7e','%15.7e','%15.7e','%15.7e','%15.7e') ,  newline='\n')



