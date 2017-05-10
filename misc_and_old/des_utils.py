import numpy as nm
import matplotlib.pyplot as plt
import healpy as hp






def plot_hist1(h1,h2,fig, **kwargs):
    """
    def plot_hist1(h1,h2,fig=figid, **kwargs):
        plt.figure(fig)
        plt.bar(h2[:-1], h1, width=h2[1]-h2[0], **kwargs)
    """

    plt.figure(fig)
    plt.bar(h2[:-1], h1, width=h2[1]-h2[0], **kwargs)

def cat2hpx_cg(ra, dec, nside):
    """
    def cat2hpx_cg(ra, dec, nside):
        R=hp.Rotator(coord='cg')
        phi, theta=R(ra, dec, lonlat=True)
        theta=90-theta
        pixs=hp.ang2pix(nside, theta*nm.pi/180, phi*nm.pi/180, nest=True)
        countmap1=nm.zeros((12*nside**2))
        for idx in pixs:
            countmap1[idx]+=1
        
        return countmap1
    """
    
    R=hp.Rotator(coord='cg')
    phi, theta=R(ra, dec, lonlat=True)
    theta=90-theta

    pixs=hp.ang2pix(nside, theta*nm.pi/180, phi*nm.pi/180, nest=True)

    countmap1=nm.zeros((12*nside**2))
    for idx in pixs:
        countmap1[idx]+=1


    return countmap1

