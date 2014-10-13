

import numpy as nm
from scipy import integrate as inte
import scipy.interpolate


class ker_shear:
    def __init__(self, cc, dNdz):
        print ' loading shear kernel'
        self.chi_star=cc.chi_de_z(1100.)
        self.cc=cc
        self.dNdz=dNdz
        
        

        self.qtab=[inte.quad (lambda chip: self.dNdz(self.cc.z_de_chi(chip))  * self.cc.dzsurdchi(chip)*   ( chip- chi   )/chip , chi, self.chi_star)[0] for chi in self.cc.chitab]

       # self.chitab=[inte.quad(lambda l: 1.0/nm.sqrt(self.Edez(l)),0.,z)[0] for z in self.ztab]


        self.spline_q_de_chi=scipy.interpolate.UnivariateSpline(self.cc.chitab, self.qtab, k=3,s=0)

    def ker(self, chi,l):

        return self.spline_q_de_chi(chi) * 3./2.* self.cc.H0**2/ self.cc.c**2 *self.cc.omegam  * chi * (1+self.cc.z_de_chi(chi))
 












#    def dNsurdz(self,z):
    
 #       z0=1.1
    
  #      if  z< z0:
   #         dndz=nm.exp(- (z-z0)**2/(2*0.8**2))
    #    else:
     #       dndz=nm.exp(- (z-z0)**2/(2*0.3**2))
      #  return dndz

#norm_N=inte.quad( lambda chi :  dNsurdz(cc.z_de_chi(chi))  * cc.dzsurdchi(chi),   0, cc.chi_de_z(1100) )[0]
#norm_Nz=inte.quad( lambda z: dNsurdz(z), 0, 1100    )[0]
