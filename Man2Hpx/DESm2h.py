

import matplotlib.pyplot as plt
import healpy as hp
import numpy as nm
import pymangle as pym
import os
from Man2Hpx import m_utils as mu




def DESm2h(tolyfile, nside_out, nside_tiny, coord='GAL',weights=['maglims'], nside_large=16, nside_ini=128, project='OPS' ):
    """
    Routine designed to convert DES mangle mask to healpix masks without having the full combined masks
    The only combined mask that is required is the "tolygon polygons" file, that represent the footprint  of the mask at the tile level
    
    Inputs:
    
    tolyfile : toylogon file
    
    Parameters:
    
    nside_out: Nside of the final Healpix mask that are produced  (typically nside_out=2048 or 4096)
    nside_tiny: Nside of the finest resolution (typically nside_tiny=16384 or 32768)
    weights: list os strings that  indicate what maps are to be produced. Typically weights=['maglims', 'bitmask', 'time', 'weight']
    nside_large: very small nside (coarse resultion to minimize time. It is faster to run polyid on a big array of points than multiple polyid calls with smaller arrays ). Typically nside_large=8 or 16
    nside_ini : smallest nside needed tobe sure to probe all the tolygon file.  One tile is .5 sq deg. that requires probing the nside=128 pixels 
    coord: coordinante system. Should be only GAL or EQU
    project: Either ACT or OPS. used to find the mangle files on deslogin

    Outputs:

    N+1 healpix maps at resolution Nside_out, in Nested format and in the desired coordinates, where N is the number of weights considered. The first map (ind =0), is the detection fraction, i.e. the fraction of tiny pixel that are not 0 inside a final pixel.
    
   
    
    """
    if coord!='GAL' and coord!='EQU':
        print  "Please indidicate a correct coordinate system. Use either GAL or EQU"
        return -1

    nmaps=len(weights)
    maps=nm.zeros((12*nside_map_final**2, nmaps+1))+hp.UNSEEN


    ####################################################################################################
    ### 1. get at least one healpix pixel center within each tile of the mask. One tile is .5 sq deg. that requires probing the nside=128 pixels 
    ######################################################


    ## load ra, dec of healpix center for nside=nside_ini and galactic coordinates
    print 'Loading healpix center coordinates for nside= %d'%nside_ini
    
    if coord=='GAL':
        
        ra,dec= mu.get_hpx_coord_GAL_inradec(nside_ini)
    if coord=='EQU':
        ra,dec= mu.get_hpx_coord_radec(nside_ini)
    
    #Load tolygon files 
    print 'loading tolygons file'
    mng=pym.Mangle(tolyfile)
    print 'done'

    ##### get the weight of the nside=128 pixels against the tolygon mask. Tell whether pixels are in the mask
    w=mng.weight(ra,dec)


    ind=nm.where(w !=0)[0]


    ### ind contains the pixels that are in the mask. to make sure that we are indeed having all the pixels that have even the tiniest overlap with the mask, we take twice the neighbouring pixels

    a=hp.get_all_neighbours(nside_ini,ind, nest=True)
    a=a[a>0]
    b=nm.concatenate([a, ind])
    b=nm.unique(b)
    c=hp.get_all_neighbours(nside_ini,b, nest=True)
    c=c[c>0]
    c=nm.concatenate([c,b])
    d=nm.unique(c)

    ### d now contains all the mid-scale pixels that we are going to consider.



    ###########################################################################################
    ### 2. consider even larger pixels to minimize io time. it is faster to run polyid on a big array of points than multiple polyid calls with smaller arrays
    ##########################################################################################
    dlarge=nm.int32(d/4**(nm.log(nside_ini/nside_large)/nm.log(2)))
    dlarge=nm.unique(dlarge)

    print "There are %d large nside=%d pixels in this mask'"%(dlarge.shape[0], nside_large)
    # dlarge now contains the ids of the nside= nside_large pixels

    ############################################################################################
    ### 3. loading the table with the mangle_run, coadd_run ad coaddtile_id
    ############################################################################################

    #### Need to add something to get the tiles.dat files


    ff=open('/Users/benoitl/Documents/Post_doc/DES/mangle_healpix/dev/y1p1_s82_tiles.dat')
    lines=ff.readlines()
    N=len(lines)-1
    dum1=[]
    for i in range(1,N+1):
        dum=lines[i].strip().split(',')
        dum1.append(dum)
    tile_tab=nm.array(dum1)
    ff.close()

    #### tile_tab contains all the info af the tiles in the mask: coaddtile_id, tilename, mangle_run, coadd_run. Those should be needed to retreive the individual tiles masks



    npix_per_fat=nm.int(4**(nm.log(nside_tiny/nside_map_final)/nm.log(2)))  ## number of tiny pixels per pixels of the final healpix map.

    kk=0
    for lpixid in dlarge:    ##### loop on the large pixels
        kk+=1
        print '%d / %d big pixels'%(kk, dlarge.shape[0])

        ## get all the pixels of the final map for the large pixel lpixid 
        fat=mu.get_sons2(nm.array([lpixid]), nside_large,nside_map_final)

        ##these are the ipix of all the tiny hpx pixels within the fat 
        sons=mu.get_sons2(fat, nside_map_final, nside_tiny)
        jmaps=nm.zeros((len(sons), nmaps+1))


        ####################################################################
        ### Need to transform back the coord of the sons into EQU to be able to match them agains the mangle mask

        Npixel=sons.shape[0]
        R=hp.Rotator(coord='gc')

        theta,phi=hp.pix2ang(nside_tiny, sons, nest=True)
        theta=nm.pi/2 -theta
        if coord=='GAL':
            ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
            ra=nm.mod(ra, 360)
        if coord=='EQU':
            ra=phi*180/nm.pi
            dec=theta*180/nm.pi

            ra=nm.mod(ra, 360)
            ra,dec=mu.perturb_radec(nside_tiny, ra,dec)



        #### Need to split the sons pixels according to which tile thy fall in.
        p=mng.polyid(ra, dec)
        a,b,c=nm.unique(p, return_index=True, return_inverse=True)

        print 'There are %d tiles in this large pixel of the sky in which we will look one by one'%(a.shape[0])
        
        #####For each tile within this large pixel lpixid interoggate the mangle masks at the position of the tiny pixels.

        for i in  range(a.shape[0]):
            t1=a[i]
            if t1!=-1: 
                ind_1= nm.where(tile_tab[:,0]=='%s'%t1)[0]  ## get the line in tile_tab that correspond to that tile

                mangle_run=tile_tab[ind_1,2][0]
                coadd_run=tile_tab[ind_1,3][0]

                # t1 is the tileid of the firsttile in which there are some pixles.
                #ind1=where(c==i) are the indices of those pixels
                indd=nm.where(c==i) 
                verbose=False

                # loading the mangle mask for that first tile. The path of the files with have to be modified when ported to the deslogin
                if verbose:
                    print 'loading file %s'%'/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/Final_masks/mangle/y1p1_s82_copy/%s/mangle/%s_holymolys_weight_i.pol'%(mangle_run, coadd_run)


                nmg1=pym.Mangle('/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/Final_masks/mangle/y1p1_s82_copy/%s/mangle/%s_holymolys_weight_i.pol'%(mangle_run, coadd_run))
                w1= nmg1.weight(ra[indd], dec[indd])
                jmaps[indd,0]=nm.float64(w1)

                #                jweight[indd]=nm.float64(w1)
                del(w1)
                it=0
                for  wei in weights:
                    it+=1
                    ### load maglim
                    nmg1.read_weights('/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/Final_masks/mangle/y1p1_s82_copy/%s/mangle/%s_molys_i.%s'%(mangle_run, coadd_run, wei))
                    w1= nmg1.weight(ra[indd], dec[indd])
                    jmaps[indd, it]=nm.float64(w1)
                    #                    jmaglim[indd]=nm.float64(w1)
                    del(w1)
                    ###fin boucle sur tile 

        jjmaps=nm.reshape(jmaps, (Npixel/npix_per_fat, npix_per_fat,nmaps+1))

        mm=nm.mean(jjmaps, axis=1, dtype=nm.float64)
        mfrac=nm.sum(jjmaps[:,:,0]>0, axis=1)/(1.0*npix_per_fat)
        maps[fat,1:]=mm[:,1:]
        maps[fat,0]=mfrac
        del( sons, fat, jjmaps,jmaps,mm, mfrac,)
        del(p, a,b,c)
    return maps











root_dir= '/Users/benoitl/Documents/Post_doc/DES/mangle_healpix/'

#man_dir='/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/coadd/20130320000001_DES0219-0541//mangle/'
#man_file='20130320000001_DES0219-0541_holymolys_weight_r.pol'

## directory of the combined mangle directory where the tolys.pol file is. 
man_dir='/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/Final_masks/mangle/y1p1_s82/'
man_file='y1p1_s82_tolys.pol'


##tolyfile='/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/Final_masks/mangle/y1p1_s82/y1p1_s82_tolys.pol'

tolyfile=man_dir+man_file

nside_large=8 ## Healpix map is done in each of these big pixels
nside_ini=128 ## Initial pixels udes to find where the mangle mask is, should be failty big, to have a dense ampling of the tolygon polygon file
nside_map_final=1024  ## Final resolution of the output healpix maps
nside_tiny=4096   ### finer scale on which the averaging is done



weights=['maglims', 'bitmask', 'time', 'weight']


#maps_G=DESm2h_gal(tolyfile, nside_map_final, nside_tiny, weights=weights, nside_large=nside_large, nside_ini=128 )
#maps_C=m2h_equ(tolyfile, nside_map_final, nside_tiny, weights=weights, nside_large=16, nside_ini=128 )









