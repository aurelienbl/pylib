

import numpy as nm
import matplotlib.pylab as plt
import scipy as si

from scipy import integrate as inte
import scipy.interpolate
import os
import spice_wr as spice
import abl_lib as abl
import pyfits as pyf
import healpy as hp
import sys


cl=abl.get_cltt(2048)

class sp_mask():
    
    def __init__(self, path, mask):
        """
        Definit une classe mask pour spice. 
        Inputs:
        path: Path de l endroit ou vont etre stockes les outputs. Doit exister avant.
        mask : masque sur le lequel on veut travailler
        
        Actions: cree les dossiers dans lesquels seront stockes des fichers, et charge aussi le masque en memoire.
       

        """

        if path[-1]!='/':
            path+='/'
        self.path=path  ### path dans lequel seront stocks tous les trucs
          ### mask 
         
        maskname=mask.split('/')[-1]
        maskname=maskname.split('.')[-2]
        self.maskname=maskname




        if os.path.isdir(self.path)==False:
            print 'You must create %s before!'%self.path
            sys.exit()
            

        if os.path.isdir(self.path+maskname)==False:
            os.mkdir(self.path+maskname)
            os.mkdir(self.path+maskname+'/temp')
            
           
        self.path_temp=self.path+'temp'
       

        if os.path.isfile(self.path+ mask.split('/')[-1])==False:
            os.system('cp %s %s'%(mask, self.path+'/'))
        self.maskf=self.path+mask.split('/')[-1]
        self.mask=hp.read_map(self.maskf)
        self.nside=hp.get_nside(self.mask)

    def show_mask(self):
        print self.mask[1]


    def plot_cormask(self, lmax, fig=None):
        """
        permet  de plotter la fonction de correlation du masque. 
        
        """

        spice.spice(self.maskf, nlmax=lmax, corfile=self.path_temp+'cor_mask.dat')
    
        cor=nm.loadtxt(self.path_temp+'cor_mask.dat')
        plt.figure(fig)
        plt.clf()
        plt.plot(cor[:,0]*180/nm.pi, cor[:,2])
        plt.title(self.maskf)





class sp_setting(sp_mask):
    def __init__(self, mask1, apodizesigma=180, thetamax=180):
        """
        definit une classe setting pour spice propre a un masque.
        """
        sp_mask.__init__(self,mask1.path, mask1.maskf  )
        
        self.apodizesigma=apodizesigma
        self.thetamax=thetamax
    


        ## make the propoer directory for this setting
        self.suff='apo%d_the%d'%(self.apodizesigma*60, self.thetamax*60)
        
        
        if os.path.isdir(self.path+self.maskname+'/'+self.suff)==False:
            os.mkdir(self.path+self.maskname+'/'+self.suff)
            
            

            os.mkdir(self.path+self.maskname+'/'+self.suff+'/kernels')
            os.mkdir(self.path+self.maskname+'/'+self.suff+'/maps')
            os.mkdir(self.path+self.maskname+'/'+self.suff+'/cls')



        self.path_ker=self.path+self.maskname+'/'+self.suff+'/kernels/'
        self.path_maps=self.path+self.maskname+'/'+self.suff+'/maps/'
        self.path_cls=self.path+self.maskname+'/'+self.suff+'/cls/'

        
        if os.path.isfile(self.path_ker+'ker.fits'):
            self.ker=pyf.open(self.path_ker+'ker.fits')[0].data[0]


        if os.path.isfile(self.path_ker+'cov.dat'):
            self.cov=nm.loadtxt(self.path_ker+'cov.dat')
            self.cov_th=nm.loadtxt(self.path_ker+'cov_th.dat')




        
    def get_kernel(self, nside,lmax):

        cltt=abl.get_cltt(2*lmax)

        map1=hp.synfast([cltt,cltt,  cltt ,cltt*0], nside,lmax=lmax, new=False, pixwin=False )
        hp.write_map(self.path_temp+'map1.fits', map1)



        spice.spice(self.path_temp+'map1.fits', nlmax=lmax, polarization=True, kernel=self.path_ker+'ker.fits', clfile=self.path_temp+'cldum.dat', apodizesigma=self.apodizesigma, thetamax=self.thetamax )
        
        self.ker=pyf.open(self.path_ker+'ker.fits')[0].data[0]
        os.system('rm '+self.path_temp+'map1.fits')
        os.system('rm '+self.path_temp+'cldum.dat')
 



    def get_cov(self, nside,lmax,N, cl):
        self.cov=nm.zeros((lmax+1, lmax+1))
        cov=self.cov.copy()

     
        for i in range(N):
            ####### simulate maps

            map1=hp.synfast(cl, nside, lmax=nside, new=False, pixwin=False )
            hp.write_map(self.path_maps+'sim_%d.fits'%i,map1)


            ##### spice the maps

        
            spice.spice(self.path_maps+'sim_%d.fits'%i, nlmax=lmax+400, clfile=self.path_cls+'cl_spice_sim_%d.dat'%i, maskfile=self.maskf, apodizesigma=self.apodizesigma, thetamax=self.thetamax)

            ######### Erase the map

            os.system('rm %s'%self.path_maps+'sim_%d.fits'%i)

        ## load the cls
        tab=nm.zeros((N, lmax+1))

        for i in range(N):
            tab[i,:]=nm.loadtxt(self.path_cls+'cl_spice_sim_%d.dat'%i)[0:lmax+1,1]

        self.tab=tab
        ## compute the cov
        self.meancl=nm.mean(tab, axis=0)
        for i in range(N):
            cov+=nm.dot(nm.matrix(tab[i,:]).T, nm.matrix(tab[i,:]))

        cov/=N
        cov-=nm.dot(nm.matrix(self.meancl).T, nm.matrix(self.meancl))
        ##cov cl, clp = E ( cl* clp) - E[cl] E[clp]
        self.cov=cov
        nm.savetxt(self.path_ker+'cov.dat', cov)
        l=nm.arange(lmax+1)
        
        cov_th=2./(2*l+1)/nm.mean(self.mask) *( cl[0:lmax+1]**2   )
        nm.savetxt(self.path_ker+'cov_th.dat', cov_th)
        




def get_cor(bincov):
    Nbins=bincov.shape[0]

    cor=bincov.copy()*0

    for i in range(Nbins):
        for j in range(Nbins):
            cor[i,j]=bincov[i,j]/nm.sqrt(bincov[i,i]*bincov[j,j])
    return cor


def plot_cor(bincor, fig=None):
    plt.figure(fig)
    plt.clf()
    plt.pcolor(    bincor - nm.diag(nm.diag(bincor))    )
    plt.colorbar()








#cl=abl.get_cltt(2048)
#mask1=sp_mask('./toto', '../mask512_gal.fits')   
#set1=sp_setting(mask1, apodizesigma=180, thetamax=180 )

#set2=sp_setting(mask1, apodizesigma=80, thetamax=80.1 )

#mask2=sp_mask('./toto', '../mask512_sv.fits') 
#set3=sp_setting(mask2, apodizesigma=15, thetamax=15.1 )


#set3.get_kernel(512, 512)
#set3.get_cov(512, 512, 100, cl)









#B=abl.create_binning_matrix(100, 512)




