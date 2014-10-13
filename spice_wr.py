import numpy as nm
import os



spice_exe='/Users/benoitl/Documents/Post_doc/Codes/PolSpice_v02-09-00/src/spice'



def spice(mapfile, nlmax, corfile= '' , clfile='', maskfile='', mapfile2='',  maskfile2='', kernel='', thetamax=180, apodizesigma=180, decouple=True, polarization=False):

    
    comm=''
    comm += ' -mapfile %s '%mapfile

    comm+= ' -nlmax %d '%nlmax

    if corfile !='':
        comm+=' -corfile %s '%corfile

    if clfile !='':
        comm+=' -clfile %s '%clfile

    if maskfile !='':
        comm+=' -maskfile %s '%maskfile


    if mapfile2 !='':
        comm+=' -mapfile2 %s '%mapfile2

    if maskfile2 !='':
        comm+=' -maskfile2 %s '%maskfile2
    
    if kernel!='':
        comm+= ' -kernelsfileout %s '%kernel


    if decouple:
        comm+= ' -decouple YES '
    
    if polarization:
        comm+=' -polarization YES'






    if thetamax!=180:
        comm+= ' -thetamax %f '%thetamax

    if apodizesigma!=180:
        comm+= ' -apodizesigma %f '%apodizesigma



    comm=spice_exe +comm
    os.system(comm)
