

import numpy as nm
from scipy import integrate as inte
import scipy.interpolate


class ker_gal:
    def __init__(self, cc, dNdz,b):
        print ' loading galaxy kernel'
        self.chi_star=cc.chi_de_z(1100.)
        self.cc=cc
        self.dNdz=dNdz
        self.b=b


    def ker(self, chi,l):
        return self.dNdz(self.cc.z_de_chi(chi))  * self.cc.dzsurdchi(chi) *self.b
 












#    def dNsurdz(self,z):
    
 #       z0=1.1
    
  #      if  z< z0:
   #         dndz=nm.exp(- (z-z0)**2/(2*0.8**2))
    #    else:
     #       dndz=nm.exp(- (z-z0)**2/(2*0.3**2))
      #  return dndz

#norm_N=inte.quad( lambda chi :  dNsurdz(cc.z_de_chi(chi))  * cc.dzsurdchi(chi),   0, cc.chi_de_z(1100) )[0]
#norm_Nz=inte.quad( lambda z: dNsurdz(z), 0, 1100    )[0]
