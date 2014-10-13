
import numpy as nm
import matplotlib.pylab as plt
import scipy as si
import pycamb as pc
from scipy import integrate as inte
import scipy.interpolate


class cambpy1:
    def __init__(self, cosmomo, NonLinear=0):

        self.omegab=cosmomo.omegab
        self.omegam=cosmomo.omegam
        self.omegal=cosmomo.omegal
        self.omegac= self.omegam-self.omegab
        self.H0=cosmomo.H0
        self.h=self.H0/100.
        self.camb_params={"omegab":self.omegab, "omegac":self.omegac, 'omegav':self.omegal, "H0":self.H0, "scalar_amp":2.1e-9, "reion__optical_depth":0.09, "reion__reionization":True, "reion__use_optical_depth":True, "reion__redshift":11, "scalar_index":0.9692, "reion__delta_redshift":1.5, "TCMB":2.726, "Num_Nu_massless":3.046, "yhe":0.24, "NonLinear":NonLinear}

        self.z_tab=[0, 0.1, 0.2, 0.3, 0.5, 0.7,1., 1.5,2,2.5,3.,3.5,4,5,7,9.,10.,15,20, 50, 70, 100, 200, 400, 600, 1000,1100,  2000]


        self.AA=pc.matter_power(redshifts=self.z_tab[::-1], maxk=1000, **self.camb_params)
        self.Pdek=scipy.interpolate.RectBivariateSpline( nm.log10(self.AA[0][:,0]*self.h),self.z_tab,(self.AA[1][:,::-1]/self.h**3), kx=3, ky=3, s=0)
#cc=cm.cosmomo(omegam=omega_m, omegal=omega_l, H0=H0)



class cambpy2:
    def __init__(self, cosmomo, NonLinear=0):

        self.omegab=cosmomo.omegab
        self.omegam=cosmomo.omegam
        self.omegal=cosmomo.omegal
        self.omegac= self.omegam-self.omegab
        self.H0=cosmomo.H0
        self.h=self.H0/100.
        self.camb_params={"omegab":self.omegab, "omegac":self.omegac, 'omegav':self.omegal, 'omegan':0.0, "H0":self.H0, "scalar_amp":2.215e-9, "reion__optical_depth":0.0925, "reion__reionization":True, "reion__use_optical_depth":True, "reion__redshift":11.37, "scalar_index":0.9619, "reion__delta_redshift":1.5, "TCMB":2.725, "Num_Nu_massless":3.046, "yhe":0.248, "NonLinear":NonLinear, "lSampleBoost":1, "AccuracyBoost":1, "lAccuracyBoost":1}

        self.z_tab=[0.001, 0.1, 0.2, 0.3,0.4, 0.5,0.6, 0.7,0.8,0.9,1., 1.1, 1.2, 1.3, 1.4, 1.5,1.7, 2,2.3, 2.5,2.7, 3.,3.5,4,5,7,9.,10.,15,20, 50, 70, 100, 200, 400, 600, 1000,1100, 2000]
        #self.z_tab=[0, 0.1, 0.2, 0.3, 0.5, 0.7, 70, 100, 200, 400, 600, 1000,1100, 2000]


        self.AA=pc.matter_power(redshifts=self.z_tab[::-1], maxk=100, logk_spacing=0.1, **self.camb_params)
        self.Pdek=scipy.interpolate.RectBivariateSpline( (self.AA[0][:,0]*self.h),self.z_tab,(self.AA[1][:,::-1]/self.h**3), kx=3, ky=3, s=0)
        
#cc=cm.cosmomo(omegam=omega_m, omegal=omega_l, H0=H0)

