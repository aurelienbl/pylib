


import numpy as nm
import matplotlib.pylab as plt
import scipy as si
#import pycamb as pc
from scipy import integrate as inte
import scipy.interpolate

class cosmomo:
    def __init__(self, omegam=0.3183, omegal=0.6817, omegak=0., omegab=0.04902,H0=67.04, sigma8=0.852, tilt=0.9619, Tcmb=2.725):
        
        self.omegam=omegam
        self.omegal=omegal
        self.omegab=omegab
        self.omegak=omegak
        self.c=299792.458   ### km/s
        self.H0=H0     # km/s/Mpc
        self.pc=3.086e13    # 1=pc=3e13 km
        self.h=self.H0/100.
        self.sigma8=sigma8
        self.tilt=tilt
        self.Tcmb=Tcmb
        self.ztab=nm.logspace(-2, nm.log10(3000), 2000)
        
        self.chitab=[inte.quad(lambda l: self.c/self.H0/nm.sqrt(self.Edez(l)),0.,z)[0] for z in self.ztab]


        self.spline_chi_de_z=scipy.interpolate.UnivariateSpline(self.ztab, self.chitab, k=3,s=0)
        self.spline_z_de_chi=scipy.interpolate.UnivariateSpline(self.chitab, self.ztab, k=3, s=0)

        self.age=1./self.H0*inte.quad(lambda l: 1./(l**2*nm.sqrt(self.Edez(1./l-1))),0,1)[0]


        
    def Edez(self,z):
        return  (self.omegak*(1.+z)**2.0+self.omegal+self.omegam*(1.0+z)**3.0)
    
    def chi_de_z(self, z):
        return self.spline_chi_de_z(z)
    def z_de_chi(self, chi):
        return self.spline_z_de_chi(chi)
  
    def dzsurdchi(self,chi):
        return  nm.sqrt(self.Edez(self.z_de_chi(chi)))/(self.c/self.H0)

    def r_de_z(self, z):
        if self.omegak>0:
            return nm.sinh(self.chi_de_z(z))
        
        if self.omegak==0:
            return self.chi_de_z(z)

        #return self.chi_de_z(z)
    def d_A(self, z):
        return self.r_de_z(z)/(1.+z)
    def d_L(self, z):
        return self.r_de_z(z)*(1.+z)

    def dVcosurdz(self,z):
        '''
        Dvcosurdz(z)=4*nm.pi *(cc.c/cc.H0)**3 * cc.chi_de_z(z)**2 /nm.sqrt(cc.Edez(z)) , en Mpc^3, seulement pour omegaak=0, pour linstant
        '''
        return 4*nm.pi *(self.c/self.H0)**3 * self.chi_de_z(z)**2 /nm.sqrt(self.Edez(z))
