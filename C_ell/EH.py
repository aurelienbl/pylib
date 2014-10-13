import numpy as nm
import matplotlib.pylab as plt
import scipy as si
from scipy import integrate as inte
import scipy.interpolate
#import cambpy as cp
import cosmomo2 as cm







class EH98:
#### Note some cosmo params are hard coded here, see EH99 for corrections


    def __init__(self, cosmomo):

        
        self.s8=0.852


        self.omegab=cosmomo.omegab
        self.omegam=cosmomo.omegam
        self.omegal=cosmomo.omegal
        self.omegac= self.omegam-self.omegab
        self.omega=self.omegam
        self.H0=cosmomo.H0
        self.h=self.H0/100.
        self.N_nu=0
        self.c=cosmomo.c

        T_cmb=2.725
	self.theta_cmb=T_cmb/2.7
	self.omhh = self.omega*self.h*self.h
        self.obhh = self.omegab*self.h*self.h


        
        self.tilt=0.9692

        self.k_equality = 0.0746*self.omhh*self.theta_cmb**(-2.)


        self.z_eq = 2.50e4*self.omhh*self.theta_cmb**(-4.) 
        z_drag = 0.313*self.omhh**(-0.419)*(1.+0.607*self.omhh**(0.674))
        z_drag = 1e0 + z_drag*self.obhh**(0.238*self.omhh**(0.223))
        z_drag = 1291e0 * self.omhh**(0.251)/  (1e0 + 0.659*self.omhh**(0.828)) * z_drag
        self.z_drag=z_drag

        self.rks=1.6*(self.obhh**0.52)*(self.omhh**0.73)*(1.+((10.4*self.omhh)**-0.95))

        self.R_drag=self.RR(self.z_drag)
        self.R_eq=self.RR(self.z_eq)

        self.k_eq=7.46e-2* self.omhh *(self.theta_cmb)**(-2)

        self.s= 2./3./self.k_eq * nm.sqrt(6./self.R_eq  ) *nm.log(  (  nm.sqrt(1+self.R_drag)+ nm.sqrt(self.R_drag+self.R_eq)  )   / (1+ nm.sqrt(self.R_eq  ))   )
        self.normm=self.sigma8(0, 8/self.h)



        self.kk=nm.logspace(-3,1, 100)

        self.z_tab=[0.001, 0.1, 0.2, 0.3,0.4, 0.5,0.6, 0.7,0.8,0.9,1., 1.1, 1.2, 1.3, 1.4, 1.5,1.7, 2,2.3, 2.5,2.7, 3.,3.5,4,5,7,9.,10.,15,20, 50, 70, 100, 200, 400, 600, 1000,1100, 2000]
        self.AA=nm.zeros((len(self.kk), len(self.z_tab)))

        for i in range(0, len(self.z_tab)):
            self.AA[:,i]=self.pdekz(self.kk, self.z_tab[i])

#        self.chitab=[inte.quad(lambda l: self.c/self.H0/nm.sqrt(self.Edez(l)),0.,z)[0] for z in self.ztab]
        self.ppdekz=scipy.interpolate.RectBivariateSpline( self.kk,self.z_tab,self.AA, kx=3, ky=3, s=0)



# self.AA=self.pc.matter_power(redshifts=self.z_tab[::-1], maxk=100, logk_spacing=0.1, **self.camb_params)
#        self.Pdek=scipy.interpolate.RectBivariateSpline( (self.AA[0][:,0]*self.h),self.z_tab,(self.AA[1][:,::-1]/self.h**3), kx=3, ky=3, s=0)
        







    def pdekz_unnorm(self,  k, z):


        ### CDM transfert
        f= 1./(1.+((k*self.s/5.4)**4)) ## eq 18
        
        q=k/13.41/self.k_eq  ## eq 10

        a1=((46.9*self.omhh)**0.670)*(1.+(32.1*self.omhh)**-0.532)
        a2=((12.0*self.omhh)**0.424)*(1.+(45.0*self.omhh)**-0.582)
        alphac=(a1**(-self.omegab/self.omega))*(a2**(-(self.omegab/self.omega)**3.))


        b1=0.944/(1.+((458.*self.omhh)**-0.708))
        b2=((0.395*self.omhh)**-0.0266)
        betac=1./(1.+b1*(((self.omegac/self.omega)**b2)-1.))

        
        C1=14.2+ 386./(1.+69.9*(q**1.08))
        C2=14.2/alphac + 386./(1.+69.9*(q**1.08))
    
        T0tilde1=  nm.log(nm.exp(1)+  1.8 * betac*q)/ ( nm.log(nm.exp(1)+  1.8*betac*q) + C1 *q**2   )  #eq 19
        T0tilde2=  nm.log(nm.exp(1)+  1.8 * betac*q)/ ( nm.log(nm.exp(1)+  1.8*betac*q) + C2 *q**2   )  #eq 19
        Tc= f* T0tilde1+ (1-f)* T0tilde2
 

        ## eq 15
        y=(1.+self.z_eq)/(1.+self.z_drag)
        g=y*(-6.*nm.sqrt(1.+y)+(2.+3.*y)*nm.log((nm.sqrt(1.+y)+1.)/(nm.sqrt(1.+y)-1.)))
        
        alphab=g*2.07*self.k_equality*self.s/((1.+self.R_drag)**0.75)

#  % Equation 23
        bn=8.41*(self.omhh**0.435)

##  % Equation 22
        ss=self.s/((1+((bn/k/self.s)**3))**(1./3.))

        ##% Equation 24
        betab=0.5+(self.omegab/self.omega) + (3.-2.*self.omegab/self.omega)*nm.sqrt(((17.2*self.omhh)**2.)+1.)

##  % Equations 19 & 21
        

        tb=nm.log(nm.exp(1)+1.8*q)/(nm.log(nm.exp(1)+1.8*q)+C1*q*q)/(1+((k*self.s/5.2)**2.))
        tb=(tb+alphab*nm.exp(-((k/self.rks)**1.4))/(1+((betab/k/self.s)**3)))*nm.sin(k*ss)/k/ss
    
 ## % Equation 8

        
        TT=self.omegab/self.omega *tb + self.omegac/self.omega *Tc













        norm=1.94e-5*self.omega**(-0.785-0.05*nm.log(self.omega))*nm.exp(-0.95*(self.tilt-1)- 0.169*(self.tilt-1)**2)
        Pdek=norm**2*(self.c /self.H0 *k)**(3+self.tilt)  * TT**2 * (2*nm.pi**2)/ k**3 

        return Pdek

    def pdekz(self, k, z):

        omegaz=self.omega*(1.+z)**3 / (   self.omegal+ self.omega*(1.+z)**3       )
        omegalz= self.omegal / (   self.omegal+ self.omega*(1.+z)**3       )

        d1z= 2.5 *omegaz/(1.+z) / ( omegaz**(4./7.)-omegalz + (1.+omegaz/2.)*(1+omegalz/70.)    )

        d0=2.5 *self.omega / ( self.omega**(4./7.)-self.omegal + (1.+self.omega/2.)*(1+self.omegal/70.)    )
        return (self.s8/self.normm)**2* self.pdekz_unnorm(k,z) * d1z**2/d0**2


    def RR(self,z):
        return 31.5 * self.obhh*self.theta_cmb**(-4)*(z/1e3)**(-1)



    def sigma8(self,z, r):

        def toto(k):
            kr=k*r
            
            res=k**3/2/nm.pi**2* self.pdekz_unnorm(k,0.0)* nm.abs(  3/kr**3*(nm.sin(kr)- kr*nm.cos(kr)))**2
            return res
       # return toto(0.002)
        return nm.sqrt(inte.quad(lambda y: toto(nm.exp(y)),  nm.log(1e-10), nm.log(1e10), limit=750   )[0])
    
    def sigma8_2(self,z, r):

        def toto(k):
            kr=k*r
            
            res=k**3/2/nm.pi**2* self.pdekz(k,0.0)* nm.abs(  3/kr**3*(nm.sin(kr)- kr*nm.cos(kr)))**2
            return res
       # return toto(0.002)
        return nm.sqrt(inte.quad(lambda y: toto(nm.exp(y)),  nm.log(1e-5), nm.log(1e2), limit=750   )[0])


#  h = cospars.h;
#  om_m = cospars.Omega_m;
#  om_b = cospars.Omega_b;
#  sigma8 = cospars.sigma_8;















class EH99:
    def __init__(self, cosmomo, omeganu):


        self.s8=cosmomo.sigma8


        self.omegab=cosmomo.omegab
        self.omegam=cosmomo.omegam
        self.omegal=cosmomo.omegal
        self.omegac= self.omegam-self.omegab
        self.H0=cosmomo.H0
        self.h=self.H0/100.
        self.N_nu=0
        self.c=cosmomo.c
   #     self.camb_params={"omegab":self.omegab, "omegac":self.omegac, 'omegav':self.omegal, "H0":self.H0, "scalar_amp":2.1e-9, "reion__optical_depth":0.09, "reion__reionization":True, "reion__use_optical_depth":True, "reion__redshift":11, "scalar_index":0.9692, "reion__delta_redshift":1.5, "TCMB":2.726, "Num_Nu_massless":3.046, "yhe":0.24, "NonLinear":NonLinear}

        self.z_tab=[0, 0.1, 0.2, 0.3, 0.5, 0.7,1., 1.5,2,2.5,3.,3.5,4,5,7,9.,10.,15,20, 50, 70, 100, 200, 400, 600, 1000,1100,  2000]

        self.T_cmb=cosmomo.Tcmb#2.725

        self.omeganu=omeganu
        self.omega=self.omegam

        self.f_nu = omeganu/self.omega
        #if self.f_nu ==0 :
         #   self.f_nu=1e-10
	self.f_baryon = self.omegab/self.omega
        

	self.theta_cmb=self.T_cmb/2.7
	self.omhh = self.omega*self.h*self.h

        self.tilt=cosmomo.tilt#0.9619


### stuff in TF_set_parameters
        self.obhh = self.omhh*self.f_baryon
# Main variables

        self.z_equality = 2.50e4*self.omhh*self.theta_cmb**(-4.) 
        self.k_equality = 0.0746*self.omhh*self.theta_cmb**(-2.)

        z_drag = 0.313*self.omhh**(-0.419)*(1.+0.607*self.omhh**(0.674))
        z_drag = 1e0 + z_drag*self.obhh**(0.238*self.omhh**(0.223))
        z_drag = 1291e0 * self.omhh**(0.251)/  (1e0 + 0.659*self.omhh**(0.828)) * z_drag
        self.z_drag=z_drag

	self.y_d = (1.+self.z_equality)/(1.+self.z_drag)

        self.R_drag = 31.5*self.obhh*self.theta_cmb**(-4.)*1000e0/(1e0 + self.z_drag)
        self.R_equality = 31.5*self.obhh*self.theta_cmb**(-4.)*1000e0/(1e0 +self.z_equality)

        self.sound_horizon = 2./3./self.k_equality*nm.sqrt(6./self.R_equality)*    nm.log(( nm.sqrt(1.+self.R_drag)+nm.sqrt(self.R_drag+self.R_equality) )  /(1.+nm.sqrt(self.R_equality)))
        self.sound_horizon2= 44.5 *nm.log(9.83/self.omhh)/nm.sqrt(1+10*self.obhh**(3./4.))
        self.p_c  = -(5.-nm.sqrt(1.+24*(1.-self.f_nu-self.f_baryon)))/4.
        self.p_cb = -(5.-nm.sqrt(1.+24*(1.-self.f_nu)))/4.
        self.f_c  = 1.-self.f_nu-self.f_baryon
        self.f_cb = 1.-self.f_nu
        self.f_nub= self.f_nu+self.f_baryon
        

        alpha_nu= (self.f_c/self.f_cb)* (2.*(self.p_c+self.p_cb)+5.)/(4.*self.p_cb+5)
        alpha_nu= alpha_nu*(1.-0.553*self.f_nub+0.126*self.f_nub**3)
        alpha_nu= alpha_nu/(1.-0.193*nm.sqrt(self.f_nu)+0.169*self.f_nu)
        alpha_nu= alpha_nu*(1.+self.y_d)**(self.p_c-self.p_cb)
        alpha_nu= alpha_nu*(1.+ (self.p_cb-self.p_c)/2. * (1.+1./(4.*self.p_c+3.)/(4.*self.p_cb+7.))/(1.+self.y_d))
        self.alpha_nu=alpha_nu
        self.beta_c=1./(1.-0.949*self.f_nub)

#### End stuff in TF_set_parameters

        self.normm=self.sigma8(0, 8/self.h)

        self.kk=nm.logspace(-3,1, 100)
        self.z_tab=[0.001, 0.1, 0.2, 0.3,0.4, 0.5,0.6, 0.7,0.8,0.9,1., 1.1, 1.2, 1.3, 1.4, 1.5,1.7, 2,2.3, 2.5,2.7, 3.,3.5,4,5,7,9.,10.,15,20, 50, 70, 100, 200, 400, 600, 1000,1100, 2000]
        self.AA=nm.zeros((len(self.kk), len(self.z_tab)))

        for i in range(0, len(self.z_tab)):
            self.AA[:,i]=self.pdekz(self.kk, self.z_tab[i])

#        self.chitab=[inte.quad(lambda l: self.c/self.H0/nm.sqrt(self.Edez(l)),0.,z)[0] for z in self.ztab]
        self.ppdekz=scipy.interpolate.RectBivariateSpline( self.kk,self.z_tab,self.AA, kx=3, ky=3, s=0)


       


    def Tmaster(self, k):
        q = k*self.theta_cmb**2/self.omhh
        gamma_eff=(nm.sqrt(self.alpha_nu) + (1.-nm.sqrt(self.alpha_nu))/	(1.+(0.43*k*self.sound_horizon2)**4))

        q_eff = q/gamma_eff
        TF_master= nm.log(nm.exp(1.0)+1.84*self.beta_c*nm.sqrt(self.alpha_nu)*q_eff)
        TF_master = TF_master/(TF_master + q_eff**2*	             (14.4 + 325./(1.+60.5*q_eff**1.11)))

      #  q_nu = 3.92*q*nm.sqrt(self.N_nu/self.f_nu)
      #  TF_master = TF_master*   (1.+(1.2*self.f_nu**(0.64)*self.N_nu**(0.3+0.6*self.f_nu))/   (q_nu**(-1.6)+q_nu**(0.8)))
        
        return TF_master


    def Dcbnu(self, k,z):
        q = k*self.theta_cmb**2/self.omhh

        #y_fs = 17.2*self.f_nu*(1.+0.488*self.f_nu**(-7./6.))*(self.N_nu*q/self.f_nu)**2
        y_fs=0
        oz=self.omega*(1.+z)**3   /(self.omegal+(1.-self.omegal-self.omega)*(1.+z)**2+self.omega*(1.+z)**3)
        olz=self.omegal    /(self.omegal+(1.-self.omegal-self.omega)*(1.+z)**2+self.omega*(1.+z)**3)
        
        D = (1.+self.z_equality)/(1.+z)*5.*oz/2.*(oz**(4./7.)-olz+(1.+oz/2.)*	 (1.+olz/70.))**(-1.)
        
        #DD0= D/((1.+self.z_equality)*5.*self.omega/2.*(self.omega**(4./7.)-self.omegal+   (1.+self.omega/2.)*(1.+self.omegal/70.))**(-1.))
        
        DD0= (1.+self.z_equality) * 5.*self.omega/2. * (   self.omega**(4./7.)- self.omegal + (1+self.omega/2.)*(1+ self.omegal/70.)    )**(-1)

        p_cb = -(5.-nm.sqrt(1.+24*(1.-self.f_nu)))/4.
    
 ##       DD_cb = (1.+(D/(1.+y_fs))**(0.7))**(-p_cb/0.7)*D**(p_cb)
        DD_cb=   ( (1+D)/(1+ y_fs) )**(-p_cb/0.7) * D**(1+p_cb)

#        DD_cbnu = ((1.-self.f_nu)**(-0.7/p_cb)  	+(D/(1.+y_fs))**(0.7))**(-p_cb/0.7)*D**(p_cb)
        DD_cbnu = (self.f_cb**(-0.7/p_cb) +(D/(1.+y_fs)))**(-p_cb/0.7)*D**(1+p_cb)
        return DD_cbnu, DD_cb, D, DD0



    
       

        
    def pdekz(self, k, z):

        return (self.s8/self.normm)**2* self.pdekz_unnorm(k,z)

    def pdekz_unnorm(self, k, z):

        norm=1.94e-5*self.omega**(-0.785-0.05*nm.log(self.omega))*nm.exp(-0.95*(self.tilt-1)- 0.169*(self.tilt-1)**2)
        #print norm**2

        #norm=self.sigma8(0, 8/self.h)

        #norm=self.s8/norm
        Tdekz=1
        D1dez=1
        D1de0=1

        growth=self.Dcbnu(k,z)
        Tcbnu= self.Tmaster(k) * growth[1]/growth[2]


        

        Tdekz=Tcbnu
        D1dez=growth[2]
        D1de0=growth[3]

        #norm =nm.sqrt(2.1e-9)
        pkz= norm**2* (self.c /self.H0 *k)**(3+self.tilt)  * self.Tmaster(k)**2  * D1dez**2/D1de0**2 *(2*nm.pi**2)/ k**3  #Tdekz**2* D1dez**2/D1de0**2   *(2*nm.pi**2)/ k**3

        return pkz





    def sigma8(self,z, r):

        def toto(k):
            kr=k*r
            
            res=k**3/2/nm.pi**2* self.pdekz_unnorm(k,0.0)* nm.abs(  3/kr**3*(nm.sin(kr)- kr*nm.cos(kr)))**2
            return res
       # return toto(0.002)
        return nm.sqrt(inte.quad(lambda y: toto(nm.exp(y)),  nm.log(1e-10), nm.log(1e10), limit=750   )[0])
    
    def sigma8_2(self,z, r):

        def toto(k):
            kr=k*r
            
            res=k**3/2/nm.pi**2* self.pdekz(k,0.0)* nm.abs(  3/kr**3*(nm.sin(kr)- kr*nm.cos(kr)))**2
            return res
       # return toto(0.002)
        return nm.sqrt(inte.quad(lambda y: toto(nm.exp(y)),  nm.log(1e-5), nm.log(1e2), limit=750   )[0])



#omega_m=0.3
#omega_l=0.7
#omega_b= 0.0445

#c=299792.458   # km/s 
#H0= 70. # km/s/Mpc

#cc=cm.cosmomo(omegam=omega_m, omegal=omega_l, omegab=omega_b, H0=H0)


##cpp=cp.cambpy2(cc)
