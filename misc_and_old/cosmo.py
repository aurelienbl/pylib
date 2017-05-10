import os
import numpy as nm
import ablcamb_lib as acl


dum_path='/Users/benoitl/Documents/Post_doc/temp/'

dir_camb_file='/Users/benoitl/Documents/Post_Doc/pylib/'#os.getcwd()+'/'#/Users/benoitl/Documents/Post_doc/Lensing/Fisher/'

par=open(dir_camb_file+'params_modele.ini')
partxt=par.read()

camb1='/Users/benoitl/Documents/Post_doc/Codes/Camb/camb_abl/camb'
camb2='/Users/benoitl/Documents/Post_doc/Codes/Camb/camb_abl/camb_abl'



class cosmo:
    def __init__(self):

       


        self.output_root='./'
        self.do_lensing='T'
        self.do_nonlinear=0
        self.use_physical='T'
      
        self.lmax_s=5000
        self.kmax_s=20000
        self.lmax_t=3000
        self.kmax_t=4000
        self.use_phys='T'
        self.ombh2= 0.022222
        self.omch2= 0.12178
        self.omnuh2=0.00062
        self.omk=0
        self.hubble=67.26
        self.w=-1
        self.pivot_s=0.002
        self.pivot_t=0.002
        self.scalar_amp=2.52221e-9
        self.scalar_spectral_index=0.9598
        self.tau=0.087

        self.lens_meth=1
        self.acc_boost=2
        self.l_acc_boost=2
        self.l_samp_boost=2


        self.N_cosmo=9
        self.thetas=self.get_theta()

        self.par_fid=[self.omnuh2, self.w, self.omk, self.omch2, self.ombh2, self.tau, self.scalar_amp, self.scalar_spectral_index, self.thetas]

    def dump_camb_params(self, params_file):
        

        newpartxt=partxt%(self.output_root,    self.do_lensing,  self.do_nonlinear, self.lmax_s,    self.kmax_s, self.lmax_t,    self.kmax_t,   self.use_physical,    self.ombh2,    self.omch2,    self.omnuh2,    self.omk,    self.hubble,    self.w,self.pivot_s,self.pivot_t, self.scalar_amp, self.scalar_spectral_index,self.tau,self.lens_meth,  self.acc_boost, self.l_acc_boost,self.l_samp_boost )
        g=open(params_file,'w')
        print>>g, newpartxt
        g.close()


#    def run_camb(params_file)



    def get_theta(self):
        old_acc_boost=self.acc_boost
        old_do_lensing=self.do_lensing
        old_l_acc_boost=self.l_acc_boost
        old_l_samp_boost=self.l_samp_boost

        old_lmax_s=self.lmax_s
        old_kmax_s=self.kmax_s

        params_file=os.getcwd()+'/'+'theta.ini'
        self.lmax_s=2
        self.kmax_s=4
        self.acc_boost=1
        self.do_lensing='F'
        self.l_acc_boost=1
        self.l_samp_boost=1

        self.dump_camb_params(params_file)
        os.system(camb1+' '+params_file+ '>truc')
        os.system("less truc|grep theta |awk {'print $5'}>truci")
        l1=nm.loadtxt('truci')

        self.lmax_s=old_lmax_s
        self.kmax_s=old_kmax_s
        self.acc_boost=old_acc_boost
        self.do_lensing=old_do_lensing
        self.l_acc_boost=old_l_acc_boost
        self.l_samp_boost=old_l_samp_boost

        
        return l1

    def get_Cls(self, lmax):
    

        self.output_root=dum_path+'temp_toto'
        self.dump_camb_params(dum_path+'/temp_camb_params.ini')
        os.system(camb1+' '+dum_path+'temp_camb_params.ini')
        UnlCls=acl.load_scalCls(dum_path+'temp_toto_scalCls.dat', lmax)
        LenCls=acl.load_lensedCls(dum_path+'temp_toto_lensedCls.dat', lmax)
        return UnlCls, LenCls


