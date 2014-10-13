import pyfits as pyf

import numpy as nm



def get_cltt(lmax):
    return fits2cl('/Users/benoitl/Documents/Post_doc/Lensing/Cls_theoriques/cltt_planck.fits', lmax)
def get_clttl(lmax):
    return fits2cl('/Users/benoitl/Documents/Post_doc/Lensing/Cls_theoriques/clttl_planck.fits', lmax)
def get_clpp(lmax):
    return fits2cl('/Users/benoitl/Documents/Post_doc/Lensing/Cls_theoriques/clpp_planck.fits', lmax)





def fits2cl(filename, lcut):
    hdulist=pyf.open(filename)
    cl=hdulist[1].data.field(0)
    return cl[0:lcut+1]




def cl2fits(cl, filename, lcut):
    """cl2fits(cl, filename, lcut)"""

    table=[pyf.Column(name='TEMPERATURE',format='1D',array=cl[0:lcut+1])]
    #print table
    tbhdu=pyf.new_table(table)
    tbhdu.writeto(filename, clobber=True)



def create_binning_matrix(Nbins, lcut):
    "create_binning_matrix(Nbins, lcut)"
### CREATION DE LA MATRICE DE BINNING ###                                       
                              
#Nbins=10                                                                       
                              

    B=nm.zeros((Nbins,lcut+1 ))
    bin1=(nm.floor(nm.linspace(2, lcut, Nbins+1)))
    tab2=nm.floor((bin1[0:Nbins]+bin1[1:Nbins+1])/2)

    tab3=nm.zeros(bin1.size, nm.int)
    for i in range(1, bin1.size):
        tab3[i]=int(bin1[i])

    tab4=nm.zeros(bin1.size-1, nm.int)
    for i in range(0, bin1.size-1):
        tab4[i]=bin1[i+1]-bin1[i]

    for i in range(0, Nbins):
        B[i,tab3[i]:tab3[i+1]]=1.0/tab4[i]

    tab4[0]=tab4[0]-2
    B[0, 0]=0
    B[0, 1]=0
    return B, bin1, tab2, tab3, tab4


def create_binning_matrix2(Nbins,lmin, lmax):
    "create_binning_matrix(Nbins, lcut)"
### CREATION DE LA MATRICE DE BINNING ###                                       
                              
#Nbins=10                                                                       
                              

    B=nm.zeros((Nbins+1,lmax+1 ))
    bin1=(nm.floor(nm.linspace(lmin, lmax, Nbins+1)))
    bin1=nm.concatenate([[0,],bin1])

    tab2=nm.floor((bin1[0:Nbins+1]+bin1[1:Nbins+2])/2)

    tab3=nm.zeros(bin1.size, nm.int)
    for i in range(1, bin1.size):
        tab3[i]=int(bin1[i])

    tab4=nm.zeros(bin1.size-1, nm.int)
    for i in range(0, bin1.size-1):
        tab4[i]=bin1[i+1]-bin1[i]

    for i in range(0, Nbins+1):
        B[i,tab3[i]:tab3[i+1]]=1.0/tab4[i]

    B[0,tab3[0]:tab3[1] ]=0

    tab4[0]=tab4[0]
    B[0, 0]=0
    B[0, 1]=0
    return B, bin1, tab2, tab3, tab4


def create_bin(tab, lmax):
    """def create_bin(tab, lmax):"""

    Nbins=tab.size-1
    tab1=nm.floor(tab)
    tab1[Nbins]=lmax
    tab1[0]=tab[0]
    tab2=nm.floor((tab1[0:Nbins]+tab1[1:Nbins+1])/2)
    tab3=nm.zeros(tab2.size, nm.int)
    for i in range(0, tab2.size):
        tab3[i]=int(tab2[i])
    return tab1,tab2, tab3, Nbins


def binnage(tab1, Nbins, tab, var):
    """def binnage(tab1, Nbins, tab, var):
    return cl_binned,var_binned"""
    cl_binned=nm.zeros(( Nbins))
    var_binned=nm.zeros(( Nbins))
    for i in range(0, Nbins):
        cl_binned[i]=nm.sum(tab[tab1[i]:tab1[i+1]])/(tab1[i+1]-tab1[i])
        var_binned[i]=nm.sum(var[tab1[i]:tab1[i+1]])/(tab1[i+1]-tab1[i])**2
    return cl_binned,var_binned

def binnage_cl(tab1, Nbins, tab ):
    cl_binned=nm.zeros(( Nbins))
    var_binned=nm.zeros(( Nbins))
    for i in range(0, Nbins):
        cl_binned[i]=nm.sum(tab[tab1[i]:tab1[i+1]])/(tab1[i+1]-tab1[i])
        var_binned[i]=nm.var(tab[tab1[i]:tab1[i+1]])/(tab1[i+1]-tab1[i])
    return cl_binned, var_binned






temp_path='/Users/benoitl/Documents/Post_doc/temp/python_temp/'
