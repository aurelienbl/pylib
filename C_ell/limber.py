
import numpy as nm
from scipy import integrate as inte

def cl_limber(l,ker1, ker2,cosmo, P ):
    return inte.quad(lambda chi:  ker1(chi,l)*ker2(chi,l) *  P((l/(chi)),cosmo.z_de_chi(chi) )/ chi**2  , 0, cosmo.chi_de_z(1100) , limit=75    )[0]

def cl_limber_zvar(l,z0, z1,ker1, ker2,cosmo, P ):
    return inte.quad(lambda chi:  ker1(chi,l)*ker2(chi,l) *  P((l/(chi)),cosmo.z_de_chi(chi) )/ chi**2  , z0, cosmo.chi_de_z(z1) , limit=75    )[0]

