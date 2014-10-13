





import numpy as nm
import matplotlib.pylab as plt
import scipy as si
from scipy import integrate as inte
import scipy.interpolate
import abl_cosmo
import abl_lib as abl
import os
import ablcamb_lib_inverted as acl

####  Program qui une matrice de Fisher.
#### Idee : pour chaque  cosmology faire 1) les derivees 2) la matrice de cov, G ou NG.




def write_matrix( A, out_file):
    nm.savetxt(out_file, A, fmt='%1.16e')


def read_matrix(in_file):
    toto=nm.loadtxt(in_file)
    return toto




class Camb_config(abl_cosmo.Cosmology):
    def __init__(self,cosmo, main_path, l_max_s_ini=7000, k_max_s_ini=10000, lmax_cosmo_der=3000,  lmax_sum=5000, lmax_cov=3000, lmax_file=6000, Nband=10, lmin=2):




        abl_cosmo.Cosmology.__init__( self,omegam=cosmo.omegam, omegal=cosmo.omegal, omegak=cosmo.omegak, omegab=cosmo.omegab,omnuh2=cosmo.omnuh2, H0=cosmo.H0, sigma8=cosmo.sigma8, tilt=cosmo.tilt, Tcmb=cosmo.Tcmb)



        self.camb1='/Users/benoitl/Documents/Post_doc/Codes/camb_april14/camb'
        self.camb2='/Users/benoitl/Documents/Post_doc/Codes/camb_april14/camb_abl'

  

        self.main_path=main_path
        if os.path.isdir(self.main_path):
            print "the directory exists youpi!!"
        else:
            print "the directory does not exist, make the dir and subdirs"
            os.mkdir(self.main_path)
            os.mkdir(self.main_path+'ngcov/')
            os.mkdir(self.main_path+'ngcov/ders')
            os.mkdir(self.main_path+'parder/')


       

        self.path_cosmo_der=self.main_path+'parder/'
        self.path_NG_der=   self.main_path+'ngcov/ders/'
        self.path_NG=self.main_path+'ngcov/'




   

    


        self.out_root_ini=self.main_path+'initial'
        self.do_lensing_ini='T'
        self.do_nonlinear_ini=0
        self.l_max_s_ini=l_max_s_ini
        self.k_max_s_ini=k_max_s_ini
        self.l_max_t_ini=3000
        self.k_max_t_ini=4000
        self.use_phys_ini='T'
        
        self.ombh2_ini= self.omegab*(self.H0/100)**2
        self.omch2_ini=   (self.omegam - self.omegab)*  (self.H0/100)**2
        self.omnuh2_ini=self.omnuh2
        self.omk_ini=0
        self.hubble_ini=self.H0
        self.w_ini=-1

        self.pivot_s_ini=0.05
        self.pivot_t_ini=0.05
        
        self.scalar_amp_ini=2.214655e-9
        self.scalar_spectral_index_ini=self.tilt
        self.tau_ini=0.09473
        
        self.lens_meth_ini=1
        self.acc_boost_ini=2
        self.l_acc_boost_ini=2
        self.l_samp_boost_ini=2


        self.eps_der=1e-4
        self.verbose=False

        self.Erase_left_right=True
    
        self.rm_temp_der=False


        self.lmax_cosmo_der=lmax_cosmo_der


        self.dum_par_file=abl.temp_path+'dum_par_file.ini'


        self.cambdic_init=self.create_camb_dic()


        self.par=open('/Users/benoitl/Documents/Post_doc/pylib/Fisher_abl/params_modele_planckbf.ini')
        self.partxt=self.par.read()
        self.partxt_ini=self.partxt


        self.l1fid=self.get_l1fid()
        



        ## paramters for bands derivstives for ng cov matrix
        
        self.lmax_sum=lmax_sum
        self.lmax_cov=lmax_cov
        self.lmax_file=lmax_file
        self.delta_l=100



        self.lmin=lmin
        self.Nband=Nband
        



        
        self.delta_p=0.02
        
        ######

       
        #### Load the clfid if they exist in the path_NG_der directory
        if os.path.isfile(self.path_NG_der+'Fiducial_scalCls.dat'):
            self.load_fid()
 







    def load_fid(self):


        self.unl_fid=acl.load_scalCls(self.path_NG_der+'Fiducial_scalCls.dat', self.lmax_file)
        self.len_fid=acl.load_lensedCls(self.path_NG_der+'Fiducial_lensedCls.dat', self.lmax_file)


    def get_l1fid(self):
        dum_par_file=abl.temp_path+'dum_par.ini'
        dic=self.create_camb_dic()
    


        l1_fid=self.get_theta(dum_par_file, dic)
        if self.verbose:
            print l1_fid
        return l1_fid




    def create_camb_dic(self):
        cambdic={}
        cambdic['output_root']=self.out_root_ini
        cambdic['do_lensing']=self.do_lensing_ini
        cambdic['do_nonlinear']=self.do_nonlinear_ini
        cambdic['lmax_s']=self.l_max_s_ini
        cambdic['kmax_s']=self.k_max_s_ini
        cambdic['lmax_t']=self.l_max_t_ini
        cambdic['kmax_t']=self.k_max_t_ini
        cambdic['use_physical']=self.use_phys_ini
        cambdic['ombh2']=self.ombh2_ini
        cambdic['omch2']=self.omch2_ini
        cambdic['omnuh2']=self.omnuh2_ini
        cambdic['omm']=self.omk_ini
        cambdic['hubble']=self.hubble_ini
        cambdic['w']=self.w_ini
        #    cambdic['omb']=omb
        #    cambdic['omcdm']=omcdm
        #    cambdic['oml']=oml
        #    cambdic['omnu']=omnu
        cambdic['pivot_s']=self.pivot_s_ini
        cambdic['pivot_t']=self.pivot_t_ini
        
        cambdic['scalar_amp']=self.scalar_amp_ini
        cambdic['scalar_spectral_index']=self.scalar_spectral_index_ini
        cambdic['tau']=self.tau_ini
        cambdic['lens_meth']=self.lens_meth_ini
        cambdic['acc_boost']=self.acc_boost_ini
        cambdic['l_acc_boost']=self.l_acc_boost_ini
        cambdic['l_samp_boost']=self.l_samp_boost_ini

        return cambdic





    def reinit_cambdic(self, cambdic):
        cambdic['output_root']=self.out_root_ini
        cambdic['do_lensing']=self.do_lensing_ini
        cambdic['do_nonlinear']=self.do_nonlinear_ini
        cambdic['lmax_s']=self.l_max_s_ini
        cambdic['kmax_s']=self.k_max_s_ini
        cambdic['lmax_t']=self.l_max_t_ini
        cambdic['kmax_t']=self.k_max_t_ini
        cambdic['use_physical']=self.use_phys_ini
        cambdic['ombh2']=self.ombh2_ini
        cambdic['omch2']=self.omch2_ini
        cambdic['omnuh2']=self.omnuh2_ini
        cambdic['omm']=self.omk_ini
        cambdic['hubble']=self.hubble_ini
        cambdic['w']=self.w_ini
        #    cambdic['omb']=self.omb
        #    cambdic['omcdm']=self.omcdm
        #    cambdic['oml']=self.oml
        #    cambdic['omnu']=self.omnu
        cambdic['pivot_s']=self.pivot_s_ini
        cambdic['pivot_t']=self.pivot_t_ini
        
        cambdic['scalar_amp']=self.scalar_amp_ini
        cambdic['scalar_spectral_index']=self.scalar_spectral_index_ini
        cambdic['tau']=self.tau_ini
        cambdic['lens_meth']=self.Lens_meth_ini
        cambdic['acc_boost']=self.acc_boost_ini
        cambdic['l_acc_boost']=self.l_acc_boost_ini
        cambdic['l_samp_boost']=self.l_samp_boost_ini
        return cambdic
        






    def dump_camb_params(self,params_file,cambdic):
   


        newpartxt=self.partxt%(cambdic['output_root'],    cambdic['do_lensing'],  cambdic['do_nonlinear'], cambdic['lmax_s'],cambdic['kmax_s'], cambdic['lmax_t'],    cambdic['kmax_t'],   cambdic['use_physical'],    cambdic['ombh2'],    cambdic['omch2'],    cambdic['omnuh2'],    cambdic['omm'],    cambdic['hubble'],    cambdic['w'],cambdic['pivot_s'],cambdic['pivot_t'], cambdic['scalar_amp'], cambdic['scalar_spectral_index'],cambdic['tau'],cambdic['lens_meth'],  cambdic['acc_boost'], cambdic['l_acc_boost'],cambdic['l_samp_boost'])
        g=open(params_file,'w')
        print>>g, newpartxt
        g.close()



        

    def get_theta(self,par_file, dic):
        dic['lmax_s']=2
        dic['kmax_s']=4
        dic['output_root']=abl.temp_path+'get_theta'

        self.dump_camb_params(par_file, dic)
        os.system(self.camb1+' '+par_file+ '>truc')
        os.system("less truc|grep \*theta |awk {'print $3'}>truci")
        l1=nm.loadtxt('truci')[0]
        os.system('rm truc truci ')
        os.system('rm %s'%abl.temp_path+'get_theta_*')
        return l1







    def get_h_forfixed_theta(self,par_ini, delta_par, dic, par_file):
    
    

        dic[par_ini]=self.cambdic_init[par_ini]+delta_par
        l1=self.get_theta(par_file, dic)
        if self.verbose:
           print l1
        old_h=self.hubble_ini
        old_old_h=50.0
        if self.verbose:
            print l1, self.l1fid

        if l1==self.l1fid:
            new_h=self.hubble_ini
        while nm.abs(l1/self.l1fid-1) >self.eps_der:

            if l1 >self.l1fid:

                new_h=old_h-nm.abs(old_h-old_old_h)/2.0

            if l1<self.l1fid:
                new_h=old_h+nm.abs(old_h-old_old_h)/2.
            if self.verbose:
                print 'new_h=',new_h
            dic['hubble']=new_h
            l1=self.get_theta(par_file, dic)
            if self.verbose:
                print l1, self.l1fid

            if l1>self.l1fid:
                old_h=new_h
            if l1<self.l1fid:
                old_old_h=old_h
                old_h=new_h


        print "#######################"
        print "new h=", new_h
        print "error on la=",nm.abs(l1/self.l1fid-1)
        print "#######################"
        return new_h



    def get_h_for_theta_der(self,delta_theta, theta_fid, dic, par_file):

        l1=self.get_theta(par_file, dic)

        old_h=self.hubble_ini
        old_old_h=50.0


        if self.verbose:
            print l1, theta_fid, delta_theta
        if l1==(theta_fid+delta_theta):
            new_h=self.hubble_ini
        while nm.abs(l1/(theta_fid+delta_theta)-1) >self.eps_der:

            if l1 >theta_fid+delta_theta:

                new_h=old_h-nm.abs(old_h-old_old_h)/2.0

            if l1<theta_fid+delta_theta:
                new_h=old_h+nm.abs(old_h-old_old_h)/2.
            if self.verbose:
                print 'new_h=',new_h
            dic['hubble']=new_h
            dum_par_file=abl.temp_path+'dum_par.ini'
            self.dump_camb_params(dum_par_file,dic)
            l1=self.get_theta(par_file, dic)


            if self.verbose:
                print l1, theta_fid, theta_fid+delta_theta

            if l1>theta_fid+delta_theta:

                old_h=new_h
            if l1<theta_fid+delta_theta:
                old_old_h=old_h
                old_h=new_h
        print "#######################"
        print "new h=", new_h
        print "error on la=",nm.abs(l1/(theta_fid+delta_theta)-1)
        print "#######################"
        return new_h











    
    
    def get_der(self,param, param_ini, delta_param, lmax):
        

        cambdic1=self.cambdic_init.copy()


        if param == 'tau':
            print 'Computing %s derivative (left side)'%param

            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_left'
            cambdic1[param]=param_ini-delta_param
            params_file_left=self.path_cosmo_der+'par_'+param+'_left.ini'
            self.dump_camb_params(params_file_left,cambdic1)
            os.system(self.camb1+' '+ params_file_left)


            print 'Computing %s derivative (right side)'%param
            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_right'
            cambdic1[param]=param_ini+delta_param 
            params_file_right=self.path_cosmo_der+'par_'+param+'_right.ini'
            self.dump_camb_params(params_file_right,cambdic1)
            os.system(self.camb1+' '+ params_file_right)



        elif param=='scalar_amp':
            print 'Computing %s derivative (left side)'%param
            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_left'
            cambdic1['scalar_amp']=param_ini-delta_param
            cambdic1['tau']=self.tau_ini-0.5*nm.log( param_ini/(param_ini-delta_param ))
            params_file_left=self.path_cosmo_der+'par_'+param+'_left.ini'
            self.dump_camb_params(params_file_left,cambdic1)
            os.system(self.camb1+' '+ params_file_left)

            
            print 'Computing %s derivative (right side)'%param

            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_right'
            cambdic1['scalar_amp']=param_ini+delta_param
            cambdic1['tau']=self.tau_ini-0.5*nm.log(param_ini /(param_ini+delta_param ))
            params_file_right=self.path_cosmo_der+'par_'+param+'_right.ini'
            self.dump_camb_params(params_file_right,cambdic1)
            os.system(self.camb1+' '+ params_file_right)


        elif param=='theta':
            print 'Computing %s derivative (left side)'%param
            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_left'
            h_left=self.get_h_for_theta_der(-delta_param, param_ini,self.cambdic_init.copy() , self.dum_par_file)
            cambdic1['hubble']=h_left
            params_file_left=self.path_cosmo_der+'par_'+param+'_left.ini'
            self.dump_camb_params(params_file_left,cambdic1)
            os.system(self.camb1+' '+ params_file_left)

            print 'Computing %s derivative (right side)'%param

            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_right'
            h_right=self.get_h_for_theta_der(delta_param, param_ini,self.cambdic_init.copy() , self.dum_par_file)
            cambdic1['hubble']=h_right
            params_file_right=self.path_cosmo_der+'par_'+param+'_right.ini'
            self.dump_camb_params(params_file_right,cambdic1)
            os.system(self.camb1+' '+ params_file_right)




        elif param =='omm':
            print 'Computing %s derivative (left side)'%param
            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_left1'
            cambdic1[param]=param_ini-delta_param
            h_left1=self.get_h_forfixed_theta(param, -delta_param,self.cambdic_init.copy(), self.dum_par_file )
            cambdic1['hubble']=h_left1
            params_file_left1=self.path_cosmo_der+'par_'+param+'_left1.ini'
            self.dump_camb_params(params_file_left1,cambdic1)
            os.system(self.camb1+' '+ params_file_left1)


            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_right1'
            cambdic1[param]=param_ini+delta_param 
            h_right1=self.get_h_forfixed_theta(param, delta_param,self.cambdic_init.copy(), self.dum_par_file )
            cambdic1['hubble']=h_right1
            params_file_right1=self.path_cosmo_der+'par_'+param+'_right1.ini'
            self.dump_camb_params(params_file_right1,cambdic1)
            os.system(self.camb1+' '+ params_file_right1)

            print 'Computing %s derivative (right side)'%param
            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_left10'
            cambdic1[param]=param_ini-delta_param/10
            h_left10=self.get_h_forfixed_theta(param, -delta_param/10,self.cambdic_init.copy(), self.dum_par_file )
            cambdic1['hubble']=h_left10
            params_file_left10=self.path_cosmo_der+'par_'+param+'_left10.ini'
            self.dump_camb_params(params_file_left10,cambdic1)
            os.system(self.camb1+' '+ params_file_left10)


            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_right10'
            cambdic1[param]=param_ini+delta_param/10 
            h_right10=self.get_h_forfixed_theta(param, delta_param/10,self.cambdic_init.copy(), self.dum_par_file )
            cambdic1['hubble']=h_right10
            params_file_right10=self.path_cosmo_der+'par_'+param+'_right10.ini'
            self.dump_camb_params(params_file_right10,cambdic1)
            os.system(self.camb1+' '+ params_file_right10)


        else:
            print 'Computing %s derivative (left side)'%param
            
            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_left'
            cambdic1[param]=param_ini-delta_param
            h_left=self.get_h_forfixed_theta(param, -delta_param,self.cambdic_init.copy(), self.dum_par_file )
            cambdic1['hubble']=h_left
            params_file_left=self.path_cosmo_der+'par_'+param+'_left.ini'
            self.dump_camb_params(params_file_left,cambdic1)
            os.system(self.camb1+' '+ params_file_left)

            print 'Computing %s derivative (right side)'%param
            cambdic1['output_root']=self.path_cosmo_der+'der_'+param+'_right'
            cambdic1[param]=param_ini+delta_param 
            h_right=self.get_h_forfixed_theta(param, delta_param,self.cambdic_init.copy(), self.dum_par_file )
            cambdic1['hubble']=h_right
            params_file_right=self.path_cosmo_der+'par_'+param+'_right.ini'
            self.dump_camb_params(params_file_right,cambdic1)
            os.system(self.camb1+' '+ params_file_right)


        if param=='omm':


            l_unl1=acl.load_scalCls(self.path_cosmo_der+'der_'+param+'_left1_scalCls.dat', lmax)
            r_unl1=acl.load_scalCls(self.path_cosmo_der+'der_'+param+'_right1_scalCls.dat', lmax)

            l_len1=acl.load_lensedCls(self.path_cosmo_der+'der_'+param+'_left1_lensedCls.dat', lmax)
            r_len1=acl.load_lensedCls(self.path_cosmo_der+'der_'+param+'_right1_lensedCls.dat', lmax)

            l_unl10=acl.load_scalCls(self.path_cosmo_der+'der_'+param+'_left10_scalCls.dat', lmax)
            r_unl10=acl.load_scalCls(self.path_cosmo_der+'der_'+param+'_right10_scalCls.dat', lmax)

            l_len10=acl.load_lensedCls(self.path_cosmo_der+'der_'+param+'_left10_lensedCls.dat', lmax)
            r_len10=acl.load_lensedCls(self.path_cosmo_der+'der_'+param+'_right10_lensedCls.dat', lmax)


            du, dl=(r_unl1-r_unl10+l_unl10-l_unl1 )/(9./5.*delta_param), (r_len1-r_len10+l_len10-l_len1 )/(9./5.*delta_param)

   

           
        else:
            l_unl=acl.load_scalCls(self.path_cosmo_der+'der_'+param+'_left_scalCls.dat', lmax)
            r_unl=acl.load_scalCls(self.path_cosmo_der+'der_'+param+'_right_scalCls.dat', lmax)

            l_len=acl.load_lensedCls(self.path_cosmo_der+'der_'+param+'_left_lensedCls.dat', lmax)
            r_len=acl.load_lensedCls(self.path_cosmo_der+'der_'+param+'_right_lensedCls.dat', lmax)


            du,dl= (r_unl-l_unl)/(2*delta_param), (r_len-l_len)/(2*delta_param)
           

        if self.Erase_left_right:

            os.system('rm %s'% self.path_cosmo_der+'der_'+param+'_left*')
            os.system('rm %s'% self.path_cosmo_der+'der_'+param+'_right*')
                
          


        return du, dl




   

    def get_all_der(self):




        a= self.get_der('omnuh2', self.omnuh2_ini, delta_omnuh2, self.lmax_cosmo_der)
       # write_matrix( a[0], self.path_cosmo_der+'d_ul_omnuh2.dat')
       # write_matrix( a[1],self.path_cosmo_der+'d_l_omnuh2.dat')
        
        b= self.get_der('w', self.w_ini, delta_w, self.lmax_cosmo_der)
       # write_matrix( b[0], self.path_cosmo_der+'d_ul_w.dat')
       # write_matrix( b[1], self.path_cosmo_der+'d_l_w.dat')
            
        c= self.get_der('omm', self.omk_ini, delta_omm, self.lmax_cosmo_der)
       # write_matrix( c[0], self.path_cosmo_der+'d_ul_omm.dat')
       # write_matrix( c[1], self.path_cosmo_der+'d_l_omm.dat')
            
        d= self.get_der('omch2', self.omch2_ini, delta_omch2, self.lmax_cosmo_der)
       # write_matrix( d[0], self.path_cosmo_der+'d_ul_omch2.dat')
       # write_matrix( d[1], self.path_cosmo_der+'d_l_omch2.dat')
            
        e= self.get_der('ombh2', self.ombh2_ini, delta_ombh2, self.lmax_cosmo_der)
       # write_matrix( e[0], self.path_cosmo_der+'d_ul_ombh2.dat')
       # write_matrix( e[1], self.path_cosmo_der+'d_l_ombh2.dat')


        f= self.get_der('tau', self.tau_ini, delta_tau, self.lmax_cosmo_der)
       # write_matrix( f[0], self.path_cosmo_der+'d_ul_tau.dat')
       # write_matrix( f[1], self.path_cosmo_der+'d_l_tau.dat')

        g= self.get_der('scalar_amp', self.scalar_amp_ini, delta_scalar_amp, self.lmax_cosmo_der)
       # write_matrix( g[0], self.path_cosmo_der+'d_ul_scalar_amp.dat')
       # write_matrix( g[1], self.path_cosmo_der+'d_l_scalar_amp.dat')

        h= self.get_der('scalar_spectral_index', self.scalar_spectral_index_ini, delta_scalar_spectral_index, self.lmax_cosmo_der)
       # write_matrix( h[0], self.path_cosmo_der+'d_ul_scalar_spectral_index.dat')
       # write_matrix( h[1], self.path_cosmo_der+'d_l_scalar_spectral_index.dat')


        i= self.get_der('theta', self.l1fid, delta_theta, self.lmax_cosmo_der)
       # write_matrix( i[0], self.path_cosmo_der+'d_ul_theta.dat')
       # write_matrix( i[1], self.path_cosmo_der+'d_l_theta.dat')


        UnlAllDer=nm.zeros((9,5,self.lmax_cosmo_der+1))
        LenAllDer=nm.zeros((9,4,self.lmax_cosmo_der+1))

        UnlAllDer[0,:,:]=a[0][:,:]
        LenAllDer[0,:,:]=a[1][:,:]

        UnlAllDer[1,:,:]=b[0][:,:]
        LenAllDer[1,:,:]=b[1][:,:]
        UnlAllDer[2,:,:]=c[0][:,:]
        LenAllDer[2,:,:]=c[1][:,:]
        UnlAllDer[3,:,:]=d[0][:,:]
        LenAllDer[3,:,:]=d[1][:,:]
        UnlAllDer[4,:,:]=e[0][:,:]
        LenAllDer[4,:,:]=e[1][:,:]
        UnlAllDer[5,:,:]=f[0][:,:]
        LenAllDer[5,:,:]=f[1][:,:]
        UnlAllDer[6,:,:]=g[0][:,:]
        LenAllDer[6,:,:]=g[1][:,:]
        UnlAllDer[7,:,:]=h[0][:,:]
        LenAllDer[7,:,:]=h[1][:,:]
        UnlAllDer[8,:,:]=i[0][:,:]
        LenAllDer[8,:,:]=i[1][:,:]

        return  UnlAllDer, LenAllDer





    def return_all_der(self ):

        if os.path.isfile(self.path_cosmo_der+'UnlAllDer.npy')==0:
#        if (os.system('ls '+path+'UnlAllDer.npy')!=0):

            UnlAllDer, LenAllDer=self.get_all_der()
            self.dUnldomnuh2=UnlAllDer[0] 
            self.dUnldw=UnlAllDer[1]
            self.dUnldomk=UnlAllDer[2]
            self.dUnldomch2=UnlAllDer[3]
            self.dUnldombh2=UnlAllDer[4]
            self.dUnldtau=UnlAllDer[5]
            self.dUnldAs=UnlAllDer[6]
            self.dUnldns=UnlAllDer[7]
            self.dUnldthetas=UnlAllDer[8]

            self.dLendomnuh2=LenAllDer[0]
            self.dLendw=LenAllDer[1]
            self.dLendomk=LenAllDer[2]
            self.dLendomch2=LenAllDer[3]
            self.dLendombh2=LenAllDer[4]
            self.dLendtau=LenAllDer[5]
            self.dLendAs=LenAllDer[6]
            self.dLendns=LenAllDer[7]
            self.dLendthetas=LenAllDer[8]
            #write_matrix(UnlAllDer, path+'UnlAllDer.dat')
            #write_matrix(LenAllDer, path+'LenAllDer.dat')
            nm.save(self.path_cosmo_der+'UnlAllDer', UnlAllDer)
            nm.save(self.path_cosmo_der+'LenAllDer', LenAllDer)
        else:
            #UnlAllDer=read_matrix(path+'UnlAllDer.dat')
            #LenAllDer=read_matrix(path+'LenAllDer.dat')
            UnlAllDer=nm.load(self.path_cosmo_der+'UnlAllDer.npy')
            LenAllDer=nm.load(self.path_cosmo_der+'LenAllDer.npy')

            self.dUnldomnuh2=UnlAllDer[0] 
            self.dUnldw=UnlAllDer[1]
            self.dUnldomk=UnlAllDer[2]
            self.dUnldomch2=UnlAllDer[3]
            self.dUnldombh2=UnlAllDer[4]
            self.dUnldtau=UnlAllDer[5]
            self.dUnldAs=UnlAllDer[6]
            self.dUnldns=UnlAllDer[7]
            self.dUnldthetas=UnlAllDer[8]

            self.dLendomnuh2=LenAllDer[0]
            self.dLendw=LenAllDer[1]
            self.dLendomk=LenAllDer[2]
            self.dLendomch2=LenAllDer[3]
            self.dLendombh2=LenAllDer[4]
            self.dLendtau=LenAllDer[5]
            self.dLendAs=LenAllDer[6]
            self.dLendns=LenAllDer[7]
            self.dLendthetas=LenAllDer[8]



















    def NG_cov_driver(self):
        
        

            
        


        ##############################
        ###  Prepare fiducial run
        ##############################
        cambdic_ini=self.create_camb_dic()
        cambdic=cambdic_ini
        params_file=self.path_NG_der+'Fiducial.ini'
        cambdic['output_root']=self.path_NG_der+'Fiducial'

        self.dump_camb_params(params_file,cambdic)
        os.system(self.camb1+' '+ params_file)

        self.load_fid()



        dum_par_file=self.path_NG_der+'dum_par.ini'



        ##lmax_sum=10*50+2


        

        #B=abl.create_binning_matrix2(self.Nband, 1,lmax_sum)
        B=abl.create_binning_matrix2(self.Nband, self.lmin,self.lmax_sum)
        
        C=B[3][1:]
        
        dl=C[1:]-C[:-1]

        

        unl_fid=acl.load_scalCls(self.path_NG_der+'Fiducial_scalCls.dat', self.lmax_file)
        unl=unl_fid.copy()


        for i in range(0, self.Nband):
            print  "NG_cov_driver, doing band %d / %d"%(i,self.Nband)
            unl=unl_fid.copy()


            unl[0][C[i]+1:C[i+1]+1]*=nm.exp(-self.delta_p)
            unl[1][C[i]+1:C[i+1]+1]*=nm.exp(-self.delta_p)
            unl[2][C[i]+1:C[i+1]+1]*=nm.exp(-self.delta_p)


            cambdic['output_root']=self.path_NG_der+'/band_%d_left'%i
            params_file_i_left=self.path_NG_der+'/params_band_%d_left.ini'%i
            self.dump_camb_params(params_file_i_left,cambdic)

            acl.save_scalCls(unl,self.path_NG_der+'/band_%d_left_scalCls.dat'%i ,self.lmax_file)


            os.system(self.camb2+' '+ params_file_i_left )

            unl=unl_fid.copy()
            unl[0][C[i]+1:C[i+1]+1]*=nm.exp(self.delta_p)
            unl[1][C[i]+1:C[i+1]+1]*=nm.exp(self.delta_p)
            unl[2][C[i]+1:C[i+1]+1]*=nm.exp(self.delta_p)


            cambdic['output_root']=self.path_NG_der+'/band_%d_right'%i
            params_file_i_right=self.path_NG_der+'/params_band_%d_right.ini'%i
            self.dump_camb_params(params_file_i_right,cambdic)

            acl.save_scalCls(unl,self.path_NG_der+'/band_%d_right_scalCls.dat'%i ,self.lmax_file)


            os.system(self.camb2+' '+ params_file_i_right )

            len_left=acl.load_lensedCls(self.path_NG_der+'/band_%d_left_abl_lensedCls.dat'%i, self.lmax_file)
            len_right=acl.load_lensedCls(self.path_NG_der+'/band_%d_right_abl_lensedCls.dat'%i, self.lmax_file)

            d_len=(len_right-len_left)/(2*self.delta_p)
            
            if self.rm_temp_der:
                os.system('rm %s'%(self.path_NG_der+'/band_%d_left_abl_lensedCls.dat'%i))
                os.system('rm %s'%(self.path_NG_der+'/band_%d_right_abl_lensedCls.dat'%i))
                os.system('rm %s'%(self.path_NG_der+'/band_%d_left_scalCls.dat'%i))
                os.system('rm %s'%(self.path_NG_der+'/band_%d_right_scalCls.dat'%i))

                os.system('rm '+params_file_i_left)
                os.system('rm '+params_file_i_right)

                os.system('rm '+ self.path_NG_der+'//band_%d_right_params.ini'%i       )
                os.system('rm '+ self.path_NG_der+'//band_%d_left_params.ini'%i       )



            nm.save(self.path_NG_der+'/dlen_over_dunl_band_%d.npy'%i, d_len)
##            write_matrix(d_len, self.path_NG_der+'/dlen_over_dunl_band_%d.dat'%i)



            ##### Clpp


            unl=unl_fid.copy()
            unl[3][C[i]+1:C[i+1]+1]*=nm.exp(-self.delta_p)
    
            cambdic['output_root']=self.path_NG_der+'/band_phi_%d_left'%i
            params_file_i_left=self.path_NG_der+'/params_band_phi_%d_left.ini'%i
            self.dump_camb_params(params_file_i_left,cambdic)
        
            acl.save_scalCls(unl,self.path_NG_der+'/band_phi_%d_left_scalCls.dat'%i ,self.lmax_file)
            os.system(self.camb2+' '+ params_file_i_left )
        

            unl=unl_fid.copy()
            unl[3][C[i]+1:C[i+1]+1]*=nm.exp(self.delta_p)
        
            cambdic['output_root']=self.path_NG_der+'/band_phi_%d_right'%i
            params_file_i_right=self.path_NG_der+'/params_band_phi_%d_right.ini'%i
            self.dump_camb_params(params_file_i_right,cambdic)
        
            acl.save_scalCls(unl,self.path_NG_der+'/band_phi_%d_right_scalCls.dat'%i ,self.lmax_file)
            os.system(self.camb2+' '+ params_file_i_right )
        
        
            len_left=acl.load_lensedCls(self.path_NG_der+'/band_phi_%d_left_abl_lensedCls.dat'%i, self.lmax_file)
            len_right=acl.load_lensedCls(self.path_NG_der+'/band_phi_%d_right_abl_lensedCls.dat'%i,self.lmax_file)
        
            d_len=(len_right-len_left)/(2*self.delta_p)
        
            if self.rm_temp_der:
                os.system('rm %s'% (self.path_NG_der+'/band_phi_%d_left_abl_lensedCls.dat'%i)   )
                os.system('rm %s'% (self.path_NG_der+'/band_phi_%d_right_abl_lensedCls.dat'%i)   )
                os.system('rm %s'% (self.path_NG_der+'/band_phi_%d_left_scalCls.dat'%i)   )
                os.system('rm %s'% (self.path_NG_der+'/band_phi_%d_right_scalCls.dat'%i)   )
                
            
                os.system('rm '+params_file_i_left)
                os.system('rm '+params_file_i_right)
                
                os.system('rm '+ self.path_NG_der+'//band_phi_%d_right_params.ini'%i       )
                os.system('rm '+ self.path_NG_der+'//band_phi_%d_left_params.ini'%i       )


#            write_matrix(d_len, self.path_NG_der+'/dlen_over_dphi_band_%d.dat'%i)
            nm.save(self.path_NG_der+'/dlen_over_dphi_band_%d.npy'%i, d_len)











    def form_G_cov(self):
        

        lll=nm.arange(0, self.lmax_cov+1)
        cov_gTTTT=2./(2*lll+1.)* self.len_fid[0][0:self.lmax_cov+1]**2
        self.cov_gTTTT=nm.diag(cov_gTTTT)
        
        cov_gTTEE=2./(2*lll+1.)* self.len_fid[2][0:self.lmax_cov+1]**2
        self.cov_gTTEE=nm.diag(cov_gTTEE)
        
        cov_gTTTE=2./(2*lll+1.)* self.len_fid[2][0:self.lmax_cov+1]* self.len_fid[0][0:self.lmax_cov+1]
        self.cov_gTTTE=nm.diag(cov_gTTTE)
        
        
        cov_gEEEE=2./(2*lll+1.)* self.len_fid[1][0:self.lmax_cov+1]**2
        self.cov_gEEEE=nm.diag(cov_gEEEE)
        
        cov_gEETE=2./(2*lll+1.)* self.len_fid[1][0:self.lmax_cov+1]*self.len_fid[2][0:self.lmax_cov+1]
        self.cov_gEETE=nm.diag(cov_gEETE)
        
        cov_gTETE=1./(2*lll+1.)* self.len_fid[2][0:self.lmax_cov+1]**2+ 1./(2*lll+1.)* self.len_fid[1][0:self.lmax_cov+1]* self.len_fid[0][0:self.lmax_cov+1]
        self.cov_gTETE=nm.diag(cov_gTETE)
        
        cov_gBBBB=2./(2*lll+1.)* self.len_fid[3][0:self.lmax_cov+1]**2
        self.cov_gBBBB=nm.diag(cov_gBBBB)
    







    def form_NG_cov(self ):
        



        

        #B=abl.create_binning_matrix2(self.Nband, 1,lmax_sum)
        B=abl.create_binning_matrix2(self.Nband, self.lmin,self.lmax_sum)
        
        C=B[3][1:]
        
        dl=C[1:]-C[:-1]
        
        
        
        if os.path.isfile(self.path_NG_der+'Fiducial_scalCls.dat')==0:
            self.NG_cov_driver()

        unl=self.unl_fid.copy()
        lenf=self.len_fid.copy()




        covtttt_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covttee_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covttbb_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covttte_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        coveeee_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        coveete_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1)) 
        coveebb_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covtete_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covtebb_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covbbbb_phi=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        
        

        covBBBB_unl=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covTTBB_unl=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covEEBB_unl=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        covTEBB_unl=nm.zeros((self.lmax_cov+1, self.lmax_cov+1))
        
        

        for i in range(0, self.Nband):
            print i
            Derivtab=nm.load(self.path_NG_der+'/dlen_over_dunl_band_%d.npy'%i)[:,0:self.lmax_cov+1]
            
            sum1=nm.sum(2./(1+2*nm.arange(C[i]+1,C[i+1]+1)))
            #    sum1= nm.sum(2.0/(1.0+2* nm.arange(2+delta_l*i, 2+delta_l*(i+1))))
            covBBBB_unl+=nm.dot(nm.matrix(Derivtab[3][0:self.lmax_cov+1]).T,nm.matrix(Derivtab[3][0:self.lmax_cov+1]))*1./dl[i]**2*sum1
            
            covTTBB_unl+=nm.dot(nm.matrix(Derivtab[0][0:self.lmax_cov+1]).T, nm.matrix(Derivtab[3][0:self.lmax_cov+1]))*1./dl[i]**2*nm.sum(2./(1.0+2* nm.arange(C[i]+1,C[i+1]+1))* (unl[2]**2/(unl[0]*unl[1]))[nm.arange(C[i]+1,C[i+1]+1)])
            
            covTEBB_unl+=nm.dot(nm.matrix(Derivtab[2][0:self.lmax_cov+1]).T,nm.matrix(Derivtab[3][0:self.lmax_cov+1]))*1./dl[i]**2*sum1
            covEEBB_unl+=nm.dot(nm.matrix(Derivtab[1][0:self.lmax_cov+1]).T, nm.matrix(Derivtab[3][0:self.lmax_cov+1]))*1./dl[i]**2*sum1
    




        for i in range(0, self.Nband):
            print i
            Derivtab_phi=nm.load(self.path_NG_der+'/dlen_over_dphi_band_%d.npy'%i)[:,0:self.lmax_cov+1]
            sum1= nm.sum(2.0/(1.0+2* nm.arange(C[i]+1,C[i+1]+1)))
            covtttt_phi+=nm.dot(nm.matrix(Derivtab_phi[0]).T,nm.matrix(Derivtab_phi[0]) *1/dl[i]**2*sum1) 
            covttee_phi+=nm.dot(nm.matrix(Derivtab_phi[0]).T,nm.matrix(Derivtab_phi[1]) *1/dl[i]**2*sum1) 
            covttte_phi+=nm.dot(nm.matrix(Derivtab_phi[0]).T,nm.matrix(Derivtab_phi[2]) *1/dl[i]**2*sum1) 
            covttbb_phi+=nm.dot(nm.matrix(Derivtab_phi[0]).T,nm.matrix(Derivtab_phi[3]) *1/dl[i]**2*sum1) 
            coveeee_phi+=nm.dot(nm.matrix(Derivtab_phi[1]).T,nm.matrix(Derivtab_phi[1]) *1/dl[i]**2*sum1) 
            coveete_phi+=nm.dot(nm.matrix(Derivtab_phi[1]).T,nm.matrix(Derivtab_phi[2]) *1/dl[i]**2*sum1)  
            coveebb_phi+=nm.dot(nm.matrix(Derivtab_phi[1]).T,nm.matrix(Derivtab_phi[3]) *1/dl[i]**2*sum1) 
           
            
            covtete_phi+=nm.dot(nm.matrix(Derivtab_phi[2]).T,nm.matrix(Derivtab_phi[2]) *1/dl[i]**2*sum1)   
            covtebb_phi+=nm.dot(nm.matrix(Derivtab_phi[2]).T,nm.matrix(Derivtab_phi[3]) *1/dl[i]**2*sum1)
            
         
            covbbbb_phi+=nm.dot(nm.matrix(Derivtab_phi[3]).T,nm.matrix(Derivtab_phi[3]) *1/dl[i]**2*sum1) 




        ## Gaussian covariances
        

        self.form_G_cov()


    




        cov_ana_TTTT= covtttt_phi + self.cov_gTTTT
        cov_ana_TTEE= covttee_phi + self.cov_gTTEE
        cov_ana_TTTE= covttte_phi + self.cov_gTTTE
        cov_ana_TTBB= covttbb_phi + covTTBB_unl
      
        cov_ana_EEEE= coveeee_phi + self.cov_gEEEE
        cov_ana_EETE= coveete_phi + self.cov_gEETE
        cov_ana_EEBB= coveebb_phi + covEEBB_unl
        cov_ana_TETE= covtete_phi + self.cov_gTETE
        cov_ana_TEBB= covtebb_phi + covTEBB_unl
        cov_ana_BBBB= covbbbb_phi + self.cov_gBBBB+covBBBB_unl
      
       
        
        cov_ana_BBBB_phi=covbbbb_phi+ self.cov_gBBBB
        cov_ana_BBBB_unl=covBBBB_unl+ self.cov_gBBBB
        

        
        nm.save(self.path_NG+'cov_ana_TTTT.npy', cov_ana_TTTT)
        nm.save(self.path_NG+'cov_ana_TTEE.npy', cov_ana_TTEE)
        nm.save(self.path_NG+'cov_ana_TTTE.npy', cov_ana_TTTE)
        nm.save(self.path_NG+'cov_ana_TTBB.npy', cov_ana_TTBB)
        nm.save(self.path_NG+'cov_ana_EEEE.npy', cov_ana_EEEE)
        nm.save(self.path_NG+'cov_ana_EETE.npy', cov_ana_EETE)
        nm.save(self.path_NG+'cov_ana_EEBB.npy', cov_ana_EEBB)
        nm.save(self.path_NG+'cov_ana_TETE.npy', cov_ana_TETE)

        nm.save(self.path_NG+'cov_ana_TEBB.npy', cov_ana_TEBB)
        nm.save(self.path_NG+'cov_ana_BBBB.npy', cov_ana_BBBB)


        nm.save(self.path_NG+'cov_ana_BBBB_phi.npy', cov_ana_BBBB_phi)
        nm.save(self.path_NG+'cov_ana_BBBB_unl.npy', cov_ana_BBBB_unl)
        




    def return_NG_cov(self):

        if os.path.isfile(self.path_NG+'cov_ana_TTTT.npy')==0:
            self.form_NG_cov()



        

        self.cov_ana_TTTT=nm.load(self.path_NG+'cov_ana_TTTT.npy')
        self.cov_ana_TTEE=nm.load(self.path_NG+'cov_ana_TTEE.npy')
        self.cov_ana_TTTE=nm.load(self.path_NG+'cov_ana_TTTE.npy')
        self.cov_ana_TTBB=nm.load(self.path_NG+'cov_ana_TTBB.npy')
        self.cov_ana_EEEE=nm.load(self.path_NG+'cov_ana_EEEE.npy')
        self.cov_ana_EETE=nm.load(self.path_NG+'cov_ana_EETE.npy')
        self.cov_ana_EEBB=nm.load(self.path_NG+'cov_ana_EEBB.npy')
        self.cov_ana_TETE=nm.load(self.path_NG+'cov_ana_TETE.npy')
        self.cov_ana_TEBB=nm.load(self.path_NG+'cov_ana_TEBB.npy')
        self.cov_ana_BBBB=nm.load(self.path_NG+'cov_ana_BBBB.npy')
        
        
        self.cov_ana_BBBB_phi=nm.load(self.path_NG+'cov_ana_BBBB_phi.npy')
        self.cov_ana_BBBB_unl=nm.load(self.path_NG+'cov_ana_BBBB_unl.npy')
          







class CMB_NG_cov():
    
    def __init__(self,lmax, dir):
        if not(os.path.isdir(dir)):
            print '%s does not exist, please do something about it'%dir
            return None

        else:
            if not(os.path.isfile(dir+'cov_ana_TTTT.npy'   )):
                print "titi"
                ### create the matrices
            else:
                print "toto"

                ### load the matrices

        








conti=False
if conti:
    cc=abl_cosmo.cosmology()
    cconf=camb_config(cc)

    lmax_cosmo_der=cconf.lmax_cosmo_der
    

    cambdic_ini=cconf.create_camb_dic()
    cambdic=cambdic_ini.copy()































class Fisher2:
    
    def __init__(self, camb_conf, lmaxF):
        self.camb_conf=camb_conf
        self.lmaxF=lmaxF


        #### Load of create derivatives





class Fisher:

    def __init__(self, path, params,lmax, lmax_phi):
        """ Fisher(path,para, lmax_cmb,lmax_phi) """
        
        ## Load or create cl's
        
        if (os.system('ls '+path+'UnlCls.dat' ) == 0) and (os.system('ls '+path+'LenCls.dat' ) == 0):
        
            self.UnlCls=nm.loadtxt(path+'UnlCls.dat')[:,0:lmax+1]
            self.LenCls=nm.loadtxt(path+'LenCls.dat')[:,0:lmax+1]
       
        else:
            self.UnlCls=params.get_Cls(lmax)[0]
       
            write_matrix(self.UnlCls, path+'UnlCls.dat')
            
            self.LenCls=params.get_Cls(lmax)[1]
             
            write_matrix(self.LenCls, path+'LenCls.dat')

         #### Load or create derivatives
            
        self.der=Derivatives(params)
        self.der.return_all_der(params,path, lmax)
        
         #### Set noise CMB + Phi

        self.noiset=nm.zeros(lmax+1)
        self.noisep=nm.zeros(lmax+1)
        self.noisephi=nm.zeros(lmax_phi+1)
        
       
        self.path=path
         ### Load noisephi if it exists

        if  (os.system('ls '+ path+'rec_noise.fits') ==0):
            self.noisephi=abl.fits2cl(path+'rec_noise.fits', lmax_phi)
#       else: 
        else:
            self.noisephi=nm.zeros(lmax_phi+1)

 #           self.noisephi=alens.rec_noiseTT(self.UnlCls[0], self.noiset, lmax_phi)
  #          abl.cl2fits(self.noisephi, path+'rec_noise.fits', lmax_phi)

       # self.Binning1=create_binning_matrix2(params.Nbins,params.le2, lmax)        
        
        
         #### Load Fisher matrices if they exist in the given path.
        
         
        self.lmaxcmb=lmax
        self.lmaxphi=lmax_phi



    def set_noise_CMB(self):

        self.noiset, self.noisep=make_multi_noise(self.noise_T, self.noise_P, self.theta_FWHM, self.lmaxcmb) 
       
    def set_noise_phiTT(self):
        self.noisephi=alens.rec_noiseTT(self.UnlCls[0], self.noiset, self.lmaxphi)
        abl.cl2fits(self.noisephi, self.path+'rec_noise.fits', self.lmaxphi)


    def form_deriv_vec_cosmo_l(self,params,B,derivtab_cosmo_l ):
        Nbins=B[3]
        D_derivtab_cosmo_l=nm.tensordot(B[0],derivtab_cosmo_l.T,1 ).T
    
        DD_derivtab_l=nm.zeros((params.N_cosmo, self.nb_cl_tot_l*Nbins))

        if self.nb_cl_tot_l==4:

            DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,0,:]
            DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,1,:]
            DD_derivtab_l[:,2*Nbins:3*Nbins]=D_derivtab_cosmo_l[:,2,:]
            DD_derivtab_l[:,3*Nbins:4*Nbins]=D_derivtab_cosmo_l[:,3,:]
        
        if self.nb_cl_tot_l==3:
            if self.isTT and self.isEE and self.isBB :

                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,0,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,1,:]
                DD_derivtab_l[:,2*Nbins:3*Nbins]=D_derivtab_cosmo_l[:,2,:]
            if self.isTT and self.isEE and self.isTE :
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,0,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,1,:]
                DD_derivtab_l[:,2*Nbins:3*Nbins]=D_derivtab_cosmo_l[:,3,:]
            if self.isTT and self.isBB and self.isTE :
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,0,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,2,:]
                DD_derivtab_l[:,2*Nbins:3*Nbins]=D_derivtab_cosmo_l[:,3,:]
            if self.isEE and self.isBB and self.isTE :
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,1,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,2,:]
                DD_derivtab_l[:,2*Nbins:3*Nbins]=D_derivtab_cosmo_l[:,3,:]
        if self.nb_cl_tot_l==2:
            if self.isTT and self.isEE:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,0,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,1,:]
            if self.isTT and self.isBB:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,0,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,2,:]
            if self.isTT and self.isTE:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,0,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,3,:]
            if self.isEE and self.isBB:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,1,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,2,:]
            if self.isEE and self.isTE:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,1,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,3,:]
            if self.isBB and self.isTE:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,2,:]
                DD_derivtab_l[:,Nbins:2*Nbins]=D_derivtab_cosmo_l[:,3,:]
        if self.nb_cl_tot_l==1:
            if self.isTT:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,0,:]
            if self.isEE:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,1,:]
            if self.isBB:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,2,:]
            if self.isTE:
                DD_derivtab_l[:,0:Nbins]=D_derivtab_cosmo_l[:,3,:]
    

        return DD_derivtab_l





    def form_deriv_vec_cosmo_ul(self,params,B,derivtab_cosmo_ul):
        Nbins=B[3]
    
        D_derivtab_cosmo_ul=nm.tensordot(B[0],derivtab_cosmo_ul.T,1 ).T

        DD_derivtab_ul=nm.zeros((params.N_cosmo, self.nb_cl_tot_ul*Nbins))

        if self.nb_cl_tot_ul==3:

            DD_derivtab_ul[:,0:Nbins]=D_derivtab_cosmo_ul[:,0,:]
            DD_derivtab_ul[:,Nbins:2*Nbins]=D_derivtab_cosmo_ul[:,1,:]
            DD_derivtab_ul[:,2*Nbins:3*Nbins]=D_derivtab_cosmo_ul[:,2,:]
      
        if self.nb_cl_tot_ul==2:
            if self.isTT and self.isEE:
                DD_derivtab_ul[:,0:Nbins]=D_derivtab_cosmo_ul[:,0,:]
                DD_derivtab_ul[:,Nbins:2*Nbins]=D_derivtab_cosmo_ul[:,1,:]
            if self.isTT and self.isTE:
                DD_derivtab_ul[:,0:Nbins]=D_derivtab_cosmo_ul[:,0,:]
                DD_derivtab_ul[:,Nbins:2*Nbins]=D_derivtab_cosmo_ul[:,2,:]
            if self.isEE and self.isTE:
                DD_derivtab_ul[:,0:Nbins]=D_derivtab_cosmo_ul[:,1,:]
                DD_derivtab_ul[:,Nbins:2*Nbins]=D_derivtab_cosmo_ul[:,2,:]
        if self.nb_cl_tot_ul==1:
            if self.isTT:
                DD_derivtab_ul[:,0:Nbins]=D_derivtab_cosmo_ul[:,0,:]
            if self.isEE:
                DD_derivtab_ul[:,0:Nbins]=D_derivtab_cosmo_ul[:,1,:]
            if self.isTE:
                DD_derivtab_ul[:,0:Nbins]=D_derivtab_cosmo_ul[:,2,:]
        return DD_derivtab_ul

            
    def form_cmb_cov_binned_len_G(self,params,lmax_cmb, B,nltt, nlpp):
        Nbins=B[3]
        print 'toto'
#### Fabrication de la matrice de covariance P a partir des formules analytiques
### faire les matrices par bloc, les binner puis les assembler
        
        
    #nltt = ist.get_cmb_noise_cl(1.0*noise_T, theta_FWHM, lmax_cmb)
    #nlpp = ist.get_cmb_noise_cl(noise_P, theta_FWHM, lmax_cmb)
        ll=nm.arange(lmax_cmb+1)
        print 'titi'
        
        covTTTT=nm.diag(2.0/((2.0*ll+1))*(self.LenCls[0][0:lmax_cmb+1]+nltt)**2)
        D_covTTTT=nm.dot(nm.dot(B[0], covTTTT), B[0].T)
        D_covTTTT=nm.diag(nm.diag(D_covTTTT))
        print 'titi2'
        covEEEE=nm.diag(2.0/((2.0*ll+1))*(self.LenCls[1][0:lmax_cmb+1]+nlpp)**2)
        D_covEEEE=nm.dot(nm.dot(B[0], covEEEE), B[0].T)
        D_covEEEE=nm.diag(nm.diag(D_covEEEE))
    
        covBBBB=nm.diag(2.0/((2.0*ll+1))*(self.LenCls[3][0:lmax_cmb+1]+nlpp)**2)
        D_covBBBB=nm.dot(nm.dot(B[0], covBBBB), B[0].T)
        D_covBBBB=nm.diag(nm.diag(D_covBBBB))
    
        covTETE=nm.diag(  1.0/((2.0*ll+1))*(  (self.LenCls[0][0:lmax_cmb+1]+nltt)*(self.LenCls[1][0:lmax_cmb+1]+nlpp) + (self.LenCls[2][0:lmax_cmb+1])**2    )    )
        D_covTETE=nm.dot(nm.dot(B[0], covTETE), B[0].T)
        D_covTETE=nm.diag(nm.diag(D_covTETE))
        
        covTTEE= nm.diag(2.0/((2.0*ll+1)) *( self.LenCls[2][0:lmax_cmb+1])**2  )
        D_covTTEE=nm.dot(nm.dot(B[0], covTTEE), B[0].T)
        D_covTTEE=nm.diag(nm.diag(D_covTTEE))
        
        covTTTE= nm.diag(2.0/((2.0*ll+1))  *  (self.LenCls[0][0:lmax_cmb+1]+nltt) * ( self.LenCls[2][0:lmax_cmb+1]) )
        D_covTTTE=nm.dot(nm.dot(B[0], covTTTE), B[0].T)
        D_covTTTE=nm.diag(nm.diag(D_covTTTE))
    
        covEETE= nm.diag(2.0/((2.0*ll+1))  *  (self.LenCls[1][0:lmax_cmb+1]+nlpp) * ( self.LenCls[2][0:lmax_cmb+1]) )
        D_covEETE=nm.dot(nm.dot(B[0], covEETE), B[0].T)
        D_covEETE=nm.diag(nm.diag(D_covEETE))
    
        D_covTTBB=nm.zeros((Nbins,Nbins))
        D_covEEBB=nm.zeros((Nbins,Nbins))
        D_covBBTE=nm.zeros((Nbins,Nbins))
    

        print 'titi3'
        cov_CMB_binned=nm.zeros((self.nb_cl_tot_l*Nbins,self.nb_cl_tot_l*Nbins) )
        print 'titi4'
        if self.nb_cl_tot_l==4:
            print 'titi5'
            cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTTTT
            cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covTTEE
            cov_CMB_binned[0:Nbins, 2*Nbins:3*Nbins]=D_covTTBB
            cov_CMB_binned[0:Nbins, 3*Nbins:4*Nbins]=D_covTTTE
        
            cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covTTEE.T
            cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covEEEE
            cov_CMB_binned[Nbins:2*Nbins, 2*Nbins:3*Nbins]=D_covEEBB
            cov_CMB_binned[Nbins:2*Nbins, 3*Nbins:4*Nbins]=D_covEETE
            
            cov_CMB_binned[2*Nbins:3*Nbins, 0:Nbins]=D_covTTBB.T
            cov_CMB_binned[2*Nbins:3*Nbins, Nbins:2*Nbins]=D_covEEBB.T
            cov_CMB_binned[2*Nbins:3*Nbins, 2*Nbins:3*Nbins]=D_covBBBB
            cov_CMB_binned[2*Nbins:3*Nbins, 3*Nbins:4*Nbins]=D_covBBTE
            

            cov_CMB_binned[3*Nbins:4*Nbins, 0:Nbins]=D_covTTTE.T
            cov_CMB_binned[3*Nbins:4*Nbins, 1*Nbins:2*Nbins]=D_covEETE.T
            cov_CMB_binned[3*Nbins:4*Nbins, 2*Nbins:3*Nbins]=D_covBBTE.T
            cov_CMB_binned[3*Nbins:4*Nbins, 3*Nbins:4*Nbins]=D_covTETE
         
        if self.nb_cl_tot_l==3:
            if self.isTT and self.isEE and self.isTE :
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTTTT
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covTTEE
                cov_CMB_binned[0:Nbins, 2*Nbins:3*Nbins]=D_covTTTE
                
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covTTEE.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covEEEE
                cov_CMB_binned[Nbins:2*Nbins, 2*Nbins:3*Nbins]=D_covEETE
                

                cov_CMB_binned[2*Nbins:3*Nbins, 0:Nbins]=D_covTTTE.T
                cov_CMB_binned[2*Nbins:3*Nbins, Nbins:2*Nbins]=D_covEETE.T
                cov_CMB_binned[2*Nbins:3*Nbins, 2*Nbins:3*Nbins]=D_covTETE
            if self.isTT and self.isEE and self.isBB :
               
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTTTT
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covTTEE
                cov_CMB_binned[0:Nbins, 2*Nbins:3*Nbins]=D_covTTBB
                
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covTTEE.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covEEEE
                cov_CMB_binned[Nbins:2*Nbins, 2*Nbins:3*Nbins]=D_covEEBB
                
                cov_CMB_binned[2*Nbins:3*Nbins, 0:Nbins]=D_covTTBB.T
                cov_CMB_binned[2*Nbins:3*Nbins, Nbins:2*Nbins]=D_covEEBB.T
                cov_CMB_binned[2*Nbins:3*Nbins, 2*Nbins:3*Nbins]=D_covBBBB

            if self.isTT and self.isBB and self.isTE  :
                
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTTTT
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covTTBB
                cov_CMB_binned[0:Nbins, 2*Nbins:3*Nbins]=D_covTTTE
                
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covTTBB.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covBBBB
                cov_CMB_binned[Nbins:2*Nbins, 2*Nbins:3*Nbins]=D_covBBTE

                cov_CMB_binned[2*Nbins:3*Nbins, 0:Nbins]=D_covTTTE.T
                cov_CMB_binned[2*Nbins:3*Nbins, Nbins:2*Nbins]=D_covBBTE.T
                cov_CMB_binned[2*Nbins:3*Nbins, 2*Nbins:3*Nbins]=D_covTETE
            if self.isEE and self.isBB and self.isTE :

                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covEEEE
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covEEBB
                cov_CMB_binned[0:Nbins, 2*Nbins:3*Nbins]=D_covEETE
                
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covEEBB.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covBBBB
                cov_CMB_binned[Nbins:2*Nbins, 2*Nbins:3*Nbins]=D_covBBTE

                cov_CMB_binned[2*Nbins:3*Nbins, 0:Nbins]=D_covEETE.T
                cov_CMB_binned[2*Nbins:3*Nbins, Nbins:2*Nbins]=D_covBBTE.T
                cov_CMB_binned[2*Nbins:3*Nbins, 2*Nbins:3*Nbins]=D_covTETE

        if self.nb_cl_tot_l==2:

            if self.isTT and self.isEE:

                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTTTT
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covTTEE
            
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covTTEE.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covEEEE

            if self.isTT and self.isBB:
                
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTTTT
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covTTBB
                
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covTTBB.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covBBBB
            
            if self.isTT and self.isTE:
                
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTTTT
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covTTTE
            
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covTTTE.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covTETE
                
            if self.isEE and self.isBB:
                
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covEEEE
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covEEBB
                
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covEEBB.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covBBBB

            if self.isEE and self.isTE:

                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covEEEE
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covETE
                
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covEETE.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covTETE

            if self.isBB and self.isTE:
                
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covBBBB
                cov_CMB_binned[0:Nbins, Nbins:2*Nbins]=D_covBBTE
                
                cov_CMB_binned[Nbins:2*Nbins, 0:Nbins]=D_covBBTE.T
                cov_CMB_binned[Nbins:2*Nbins, Nbins:2*Nbins]=D_covTETE

        if self.nb_cl_tot_l==1:
            if self.isTT:
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTTTT
            if self.isEE:
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covEEEE
            if self.isBB:
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covBBBB
            if self.isTE:
                cov_CMB_binned[0:Nbins, 0:Nbins]=D_covTETE
        print 'toto7'
        return cov_CMB_binned





    def form_cmb_cov_binned_unl(self,params,lmax_cmb, B,nltt, nlpp):
        Nbins=B[3]
        ll=nm.arange(lmax_cmb+1)
        
        covTTTT_ul=nm.diag(2.0/((2.0*ll+1))*(self.UnlCls[0][0:lmax_cmb+1]+nltt)**2)
        print covTTTT_ul.shape
        print B[0].shape
        D_covTTTT_ul=nm.dot(nm.dot(B[0], covTTTT_ul), B[0].T)
        D_covTTTT_ul=nm.diag(nm.diag(D_covTTTT_ul))
        
        covEEEE_ul=nm.diag(2.0/((2.0*ll+1))*(self.UnlCls[1][0:lmax_cmb+1]+nlpp)**2)
        D_covEEEE_ul=nm.dot(nm.dot(B[0], covEEEE_ul), B[0].T)
        D_covEEEE_ul=nm.diag(nm.diag(D_covEEEE_ul))        
        
    
        covTETE_ul=nm.diag(  1.0/((2.0*ll+1))*(  (self.UnlCls[0][0:lmax_cmb+1]+nltt)*(self.UnlCls[1][0:lmax_cmb+1]+nlpp) + (self.UnlCls[2][0:lmax_cmb+1])**2    )    )
        D_covTETE_ul=nm.dot(nm.dot(B[0], covTETE_ul), B[0].T)
        D_covTETE_ul=nm.diag(nm.diag(D_covTETE_ul))
    
    
        covTTEE_ul= nm.diag(2.0/((2.0*ll+1)) *( self.UnlCls[2][0:lmax_cmb+1])**2  )
        D_covTTEE_ul=nm.dot(nm.dot(B[0], covTTEE_ul), B[0].T)
        D_covTTEE_ul=nm.diag(nm.diag(D_covTTEE_ul))
    
    
        covTTTE_ul= nm.diag(2.0/((2.0*ll+1))  *  (self.UnlCls[0][0:lmax_cmb+1]+nltt) * ( self.UnlCls[2][0:lmax_cmb+1]) )
        D_covTTTE_ul=nm.dot(nm.dot(B[0], covTTTE_ul), B[0].T)
        D_covTTTE_ul=nm.diag(nm.diag(D_covTTTE_ul))
    
    
        covEETE_ul= nm.diag(2.0/((2.0*ll+1))  *  (self.UnlCls[1][0:lmax_cmb+1]+nlpp) * ( self.UnlCls[2][0:lmax_cmb+1]) )
        D_covEETE_ul=nm.dot(nm.dot(B[0], covEETE_ul), B[0].T)
        D_covEETE_ul=nm.diag(nm.diag(D_covEETE_ul))
        
        cov_CMB_binned_ul=nm.zeros((self.nb_cl_tot_ul*Nbins,self.nb_cl_tot_ul*Nbins) )

        if self.nb_cl_tot_ul==3:

            cov_CMB_binned_ul[0:Nbins, 0:Nbins]=D_covTTTT_ul
            cov_CMB_binned_ul[0:Nbins, Nbins:2*Nbins]=D_covTTEE_ul
            cov_CMB_binned_ul[0:Nbins, 2*Nbins:3*Nbins]=D_covTTTE_ul
            
            cov_CMB_binned_ul[Nbins:2*Nbins, 0:Nbins]=D_covTTEE_ul.T
            cov_CMB_binned_ul[Nbins:2*Nbins, Nbins:2*Nbins]=D_covEEEE_ul
            cov_CMB_binned_ul[Nbins:2*Nbins, 2*Nbins:3*Nbins]=D_covEETE_ul
            
            cov_CMB_binned_ul[2*Nbins:3*Nbins, 0:Nbins]=D_covTTTE_ul.T
            cov_CMB_binned_ul[2*Nbins:3*Nbins, Nbins:2*Nbins]=D_covEETE_ul.T
            cov_CMB_binned_ul[2*Nbins:3*Nbins, 2*Nbins:3*Nbins]=D_covTETE_ul
            
        if self.nb_cl_tot_ul==2:

        
            if self.isTT and self.isEE :

                cov_CMB_binned_ul[0:Nbins, 0:Nbins]=D_covTTTT_ul
                cov_CMB_binned_ul[0:Nbins, Nbins:2*Nbins]=D_covTTEE_ul
                
                cov_CMB_binned_ul[Nbins:2*Nbins, 0:Nbins]=D_covTTEE_ul.T
                cov_CMB_binned_ul[Nbins:2*Nbins, Nbins:2*Nbins]=D_covEEEE_ul

            if self.isTT and self.isTE:
 
                cov_CMB_binned_ul[0:Nbins, 0:Nbins]=D_covTTTT_ul
                cov_CMB_binned_ul[0:Nbins, Nbins:2*Nbins]=D_covTTTE_ul
                
                cov_CMB_binned_ul[Nbins:2*Nbins, 0:Nbins]=D_covTTTE_ul.T
                cov_CMB_binned_ul[Nbins:2*Nbins, Nbins:2*Nbins]=D_covTETE_ul
                
            if self.isEE and self.isTE:
                
                cov_CMB_binned_ul[0:Nbins, 0:Nbins]=D_covEEEE_ul
                cov_CMB_binned_ul[0:Nbins, Nbins:2*Nbins]=D_covEETE_ul
                
                cov_CMB_binned_ul[Nbins:2*Nbins, 0:Nbins]=D_covEETE_ul.T
                cov_CMB_binned_ul[Nbins:2*Nbins, Nbins:2*Nbins]=D_covTETE_ul
                
        if self.nb_cl_tot_ul==1: 
            
            if self.isTT:
                cov_CMB_binned_ul[0:Nbins, 0:Nbins]=D_covTTTT_ul
            if self.isEE:
                cov_CMB_binned_ul[0:Nbins, 0:Nbins]=D_covEEEE_ul
            if self.isTE:
                cov_CMB_binned_ul[0:Nbins, 0:Nbins]=D_covTETE_ul


        return cov_CMB_binned_ul



    def form_fisher_phi_unbinned(self,params, lmax_phi, noise_phi):

        Fisher_phi_cosmo=nm.zeros((params.N_cosmo, params.N_cosmo))
    ## derivative of clphiphi with respect to cosmo params

        deriv_cosmo=nm.zeros((params.N_cosmo, lmax_phi+1))
        deriv_cosmo[0,:]=self.der.dUnldomnuh2[3,0:lmax_phi+1]
        deriv_cosmo[1,:]=self.der.dUnldw[3,0:lmax_phi+1]
        deriv_cosmo[2,:]=self.der.dUnldomk[3,0:lmax_phi+1]
        deriv_cosmo[3,:]=self.der.dUnldomch2[3,0:lmax_phi+1]
        deriv_cosmo[4,:]=self.der.dUnldombh2[3,0:lmax_phi+1]
        deriv_cosmo[5,:]=self.der.dUnldtau[3,0:lmax_phi+1]
        deriv_cosmo[6,:]=self.der.dUnldAs[3,0:lmax_phi+1]
        deriv_cosmo[7,:]=self.der.dUnldns[3,0:lmax_phi+1]
        deriv_cosmo[8,:]=self.der.dUnldthetas[3,0:lmax_phi+1]
        
        for c in range(0,params.N_cosmo):
            for cp in range(0,params.N_cosmo):
                for l in range(2,lmax_phi+1):
                    Fisher_phi_cosmo[c,cp]=Fisher_phi_cosmo[c,cp]+deriv_cosmo[c,l]*deriv_cosmo[cp,l]/(self.UnlCls[3][l]+noise_phi[l])**2*(2.0*l+1)/2.
            

        return Fisher_phi_cosmo


     





    def get_fisher(self,params):
        lmax_cmb=self.lmaxcmb
        
        self.derivtab_cosmo_ul=nm.zeros((params.N_cosmo,5,lmax_cmb+1))
        
        self.derivtab_cosmo_ul[0, :,:]=self.der.dUnldomnuh2[:, 0:lmax_cmb+1]
        
        self.derivtab_cosmo_ul[1, :,:]=self.der.dUnldw[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_ul[2, :,:]=self.der.dUnldomk[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_ul[3, :,:]=self.der.dUnldomch2[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_ul[4, :,:]=self.der.dUnldombh2[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_ul[5, :,:]=self.der.dUnldtau[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_ul[6, :,:]=self.der.dUnldAs[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_ul[7, :,:]=self.der.dUnldns[:, 0:lmax_cmb+1]
        
        self.derivtab_cosmo_ul[8, :,:]=self.der.dUnldthetas[:, 0:lmax_cmb+1]


        self.derivtab_cosmo_l=nm.zeros((params.N_cosmo,4,lmax_cmb+1))
        self.derivtab_cosmo_l[0, :,:]=self.der.dLendomnuh2[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_l[1, :,:]=self.der.dLendw[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_l[2, :,:]=self.der.dLendomk[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_l[3, :,:]=self.der.dLendomch2[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_l[4, :,:]=self.der.dLendombh2[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_l[5, :,:]=self.der.dLendtau[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_l[6, :,:]=self.der.dLendAs[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_l[7, :,:]=self.der.dLendns[:, 0:lmax_cmb+1]
        self.derivtab_cosmo_l[8, :,:]=self.der.dLendthetas[:, 0:lmax_cmb+1]

        

        self.D_derivtab_cosmo_l=nm.tensordot(self.Binning1[0],self.derivtab_cosmo_l.T,1 ).T
        self.D_derivtab_cosmo_ul=nm.tensordot(self.Binning1[0],self.derivtab_cosmo_ul.T,1 ).T
        
        self.DD_derivtab_cosmo_l=self.form_deriv_vec_cosmo_l(params,self.Binning1,self.derivtab_cosmo_l)
        self.DD_derivtab_cosmo_ul=self.form_deriv_vec_cosmo_ul(params,self.Binning1,self.derivtab_cosmo_ul )

        Cov_cmb_binned_unl=self.form_cmb_cov_binned_unl(params,lmax_cmb, self.Binning1,self.noiset, self.noisep)
        Cov_cmb_binned_len_G=self.form_cmb_cov_binned_len_G(params,lmax_cmb, self.Binning1,self.noiset, self.noisep)

        self.Fisher_cmb_cosmo_binned_unl=nm.dot(nm.dot(self.DD_derivtab_cosmo_ul,nm.array(nm.matrix(Cov_cmb_binned_unl).I)), self.DD_derivtab_cosmo_ul.T)
        self.Fisher_cmb_cosmo_binned_lenG=nm.dot(nm.dot(self.DD_derivtab_cosmo_l,nm.array(nm.matrix(Cov_cmb_binned_len_G).I)), self.DD_derivtab_cosmo_l.T)
        self.Fisher_phi_cosmo_unbinned=self.form_fisher_phi_unbinned( params,self.lmaxphi, self.noisephi)




#path_cosmo_der='Users/benoitl/Documents/Post_doc/Lensing/Fisher/'





   
delta_omm=0.02  # step for omk
delta_ombh2=0.005
delta_omch2=0.005
delta_omnuh2=0.0002   # step for omnuh2
delta_w=0.05    # step for w
delta_tau=   0.01
delta_scalar_amp= 0.02e-9
delta_scalar_spectral_index=0.01
delta_theta=0.002

