import os
import numpy as nm
import abl_lib as abl
import ablcamb_lib as acl
import pipe_tools_mbp as ptm






def rec_noiseTT(cltt, noiset_debeamed ,lmaxphi):

## take numpy arrays as input and output, even if it writes stuff on disk
    
##cltot=abl.fits2cl(clttl_th_file, lmax)+ (noise *nm.pi/180./60.)**2 /  beam**2

##cl_lens_BN=abl.fits2cl(clttl_th_file, lmax)*beam**2+(noise *nm.pi/180./60.)**2
    cltot= cltt+  noiset_debeamed
    cltot[0:2]=0
    f1=1.0/cltot
    f1[0:2]=0
    
    f2=cltt/cltot
    f2[0:2]=0
    
    f1_file=ptm.path_tmp+'f1_file.fits'
    f2_file=ptm.path_tmp+'f2_file.fits'

    abl.cl2fits(f1, f1_file, lmaxphi)
    abl.cl2fits(f2, f2_file, lmaxphi)

    cltot_file=ptm.path_tmp+'cltot.fits'
    cltt_file=ptm.path_tmp+'cltt.fits' 
    abl.cl2fits(cltot, cltot_file, lmaxphi)
   
    abl.cl2fits(cltt,cltt_file, lmaxphi)


    N0Ef=ptm.path_tmp+'N0E.fits'
    N0Bf=ptm.path_tmp+'N0B.fits'
    Af=ptm.path_tmp+'A.fits'
    ptm.launch_estim(compute_biases_only='.true.', normalize='.true.', map_file1='', nside=1024, lmax=lmaxphi,f1_file=f1_file, f2_file=f2_file, underlying_cl_file=cltt_file  ,cl_tot_file= cltot_file ,theta_FWHM=5,beam_file='',cl_estim_grad=''  ,cl_estim_curl= ''  ,def_field_grad= ''  ,def_field_curl= '',N0E= '!'+N0Ef   ,N0B='!'+N0Bf ,A= '!'+Af  ,cl_unlens=cltt_file ,namelist_filename=ptm.path_nml+'estim_biais.nml', extra_depth=[] )
    
    N0E=abl.fits2cl(N0Ef, lmaxphi)
    l=nm.arange(lmaxphi+1)
    N0E/=(l*(l+1))
    N0E[0:2]=0
    return N0E
    
    


    
#launch_estim(compute_biases_only=, normalize=, map_file1=,map_file2=, nside=, lmax=,f1_file=, f2_file=, underlying_cl_file=,cl_tot_file=,theta_FWHM=,beam_file=,cl_estim_grad=,cl_estim_curl=,def_field_grad=,def_field_curl=,N0E=,N0B=,A=,cl_unlens=,N1E=,N1B=,cl_phi_file=,N2= ,namelist_filename=, extra_depth=[] )


