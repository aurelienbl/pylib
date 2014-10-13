import numpy as nm
import pymangle as pym
import healpy as hp


def get_hpx_coord_radec(nside):
    """outputs the coordinates of the center of the healpix pixels at a given nside. coordinates are in  RA,DEC
    """

    theta, phi=hp.pix2ang(nside, nm.arange(12*nside**2), nest=True)
    theta=nm.pi/2 -theta

    ra=phi*180/nm.pi
    dec=theta*180/nm.pi

    return ra, dec



def get_hpx_coord_GAL_inradec(nside):
    """outputs the coordinates of the center of the healpix pixels at a given nside. coordinates are rotate  in  RA,DEC
    """


    theta, phi=hp.pix2ang(nside, nm.arange(12*nside**2), nest=True)
    theta=nm.pi/2 -theta

    R=hp.Rotator(coord='gc')

    ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
    ra=nm.mod(ra, 360)
    return ra, dec





def perturb_radec(nside, ra,dec):
    """ An issue with Mangle is that when a ra,dec is exactly on top of a mangle pixelization line, then polyid fails. This typically happens for healpix pixels in EQU as the Mangle simple pixelization basically lies on top of some of the healpix pixel centers. This function thus slightly perturbs the given radec (typically those of the hpx center). add a random value to the ra, generated from a Gaussian with sigma = size of the hpxpixel/ 100.
    """
    
    lenn=hp.nside2resol(nside, arcmin=True)/60  /100.
    N=len(ra)

    ra1=ra+lenn*nm.random.rand(N)
    ra1=nm.mod(ra1, 360)
   
    return ra1,dec


def get_sons2(ifat, Ni, Nf):
    tab1=nm.array(ifat, dtype=nm.int)
    ng_gen=nm.log(Nf/Ni)/nm.log(2)
    nn=4**nm.int32(ng_gen)
    tab=nm.zeros(len(tab1)*4**nm.int32(ng_gen), dtype=nm.int)

    for i in range(len(tab1)):
        tab[i*nn:i*nn+nn]=range(nn*tab1[i], nn*tab1[i]+nn)
    return tab




def get_hpx_coord_radec(nside):
    """outputs the coordinates of the center of the healpix pixels at a given nside. coordinates are in  RA,DEC
    """

    theta, phi=hp.pix2ang(nside, nm.arange(12*nside**2), nest=True)
    theta=nm.pi/2 -theta

    ra=phi*180/nm.pi
    dec=theta*180/nm.pi

    return ra, dec




def get_hpx_coord_GAL_inradec(nside):
    """outputs the coordinates of the center of the healpix pixels at a given nside. coordinates are rotate  in  RA,DEC
    """


    theta, phi=hp.pix2ang(nside, nm.arange(12*nside**2), nest=True)
    theta=nm.pi/2 -theta

    R=hp.Rotator(coord='gc')

    ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
    ra=nm.mod(ra, 360)
    return ra, dec




def perturb_radec(nside, ra,dec):
    """ An issue with Mangle is that when a ra,dec is exactly on top of a mangle pixelization line, then polyid fails. This typically happens for healpix pixels in EQU as the Mangle simple pixelization basically lies on top of some of the healpix pixel centers. This function thus slightly perturbs the given radec (typically those of the hpx center). add a random value to the ra, generated from a Gaussian with sigma = size of the hpxpixel/ 100.
    """
    
    lenn=hp.nside2resol(nside, arcmin=True)/60  /100.
    N=len(ra)



  
    ra1=ra+lenn*nm.random.rand(N)
    ra1=nm.mod(ra1, 360)

    #dec1=dec+(lenn* (nm.random.rand(N)+1))*sig2

  
    return ra1,dec



def make_mask_EQU( mangle_file, nside_out):
    """ Convert a mangle polygon file into a EQU healpix map
    """


    #### 1st step:  get the coords of the center of the hpx pixels  expressed in ra,dec

    ra,dec= get_hpx_coord_radec(nside_out)
    ra,dec=perturb_radec(nside_out, ra,dec)
    #### 2 nd step : load the mangle file with pymangle
    mng=pym.Mangle(mangle_file)

    #### 3rd step run polyid
    p=mng.polyid_and_weight(ra,dec)

    #### 4th step reproject into hpx pixel
    ind_unseen=nm.where(p[0]==-1)
    ind_seen=nm.where(p[0]!=-1)
    mask_out=nm.zeros(12*nside_out**2, dtype=nm.float128)

    mask_out[ind_unseen[0]]=hp.UNSEEN
    mask_out[ind_seen[0]]=p[1][ind_seen]
    
    return mask_out

    
    




def make_mask_GAL( mangle_file, nside_out):
    """ Convert a mangle polygon file into a GAL healpix map
    """

    #### 1st step:  get the coords of the center of the hpx GAL pixels rotatec in   expressed in ra,dec

    ra,dec= get_hpx_coord_GAL_inradec(nside_out)
    #ra,dec=perturb_radec(nside_out, ra,dec)
    #### 2 nd step : l-oad the mangle file with pymangle
    mng=pym.Mangle(mangle_file)

    #### 3rd step run polyid
    p=mng.polyid_and_weight(ra,dec)

    #### 4th step reproject into hpx pixel
    ind_unseen=nm.where(p[0]==-1)
    ind_seen=nm.where(p[0]!=-1)
    mask_out=nm.zeros(12*nside_out**2, dtype=nm.float128)

    mask_out[ind_unseen[0]]=hp.UNSEEN
    mask_out[ind_seen[0]]=p[1][ind_seen]
    
    return mask_out

    
    
