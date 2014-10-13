
import numpy as nm
import os

import matplotlib.pyplot as plt






def ninebynine(mat,par, hold=False):
    par_vec=[par.omnuh2, par.w,par.omk, par.omch2*100, par.ombh2*100, par.tau*100, par.scalar_amp,par.scalar_spectral_index*100, par.thetas*100 ]
    par_lab=['$\Omega_{nu}h^2$','$w$', '$\Omega_K$', '$\Omega_{c}h^2$','$ \Omega_{b}h^2$',' $\\tau$', '$A_s$', ' n_s$', '$100\Theta_s$' ]
    plt.figure(1)
    if ( not(hold)):
        plt.clf()
    for j in range(0,9):
        print j
        for i in range(j+1, 9):
            g=get_ellipse(par.par_fid,mat,i,j,[] )
            print (j-1+(i)*8)
            plt.subplot(8,8,j*8+i) 
            plt.subplots_adjust(wspace=0.05, hspace=0.05)
            plt.plot(g[0], g[1])
            if i==j+1:
                #plt.plot(g[0],g[1])
                plt.xlabel(par_lab[i])
                plt.ylabel(par_lab[j])
            else:
                #plt.plot(g[0],g[1])
                plt.setp(plt.gca(), xticklabels=[], yticklabels=[])
                



def sixbysix(mat,par, hold=False):
    par_vec=[par.omnuh2, par.w,par.omk, par.omch2*100, par.ombh2*100, par.tau*100, par.scalar_amp,par.scalar_spectral_index*100, par.thetas*100 ]
    par_lab=['$\Omega_{nu}h^2$','$w$', '$\Omega_K$', '$\Omega_{c}h^2$','$ \Omega_{b}h^2$','$\\tau $', '$A_s$', '$n_s$', '$100 \Theta_s$' ]
    mat1=mat.copy()
    plt.figure(1, figsize=(15,15))

    if ( not(hold)):
        plt.clf()
    for i in range(0,6):
        print i
        for j in range(i+1,6):
            g=get_ellipse(par.par_fid,mat,i+3,j+3,[0,1,2] )
            plt.subplot(5,5,i+1+(j-1)*5) 
            plt.subplots_adjust(wspace=0.1, hspace=0.1)
            plt.plot(g[0], g[1])
            
            deltax=nm.abs(nm.max(g[0])-nm.min(g[0]))/nm.min(g[0])
            deltay=nm.abs(nm.max(g[1])-nm.min(g[1]))/nm.min(g[1])
            print deltax
            plt.xticks(nm.arange(nm.min(g[0])-deltax*nm.min(g[0])*.01,nm.max(g[0])+deltax*nm.min(g[0])*0, deltax*nm.min(g[0])/3), rotation=90)
            plt.yticks(nm.arange(nm.min(g[1])-deltay*nm.min(g[1])*.01,nm.max(g[1])+deltay*nm.min(g[1])*0, deltay*nm.min(g[1])/3))
            
            if i==0 :
                plt.ylabel(par_lab[j+3])
            if j==5:
                plt.xlabel(par_lab[i+3])
          
            if i>0 and j<5:
                plt.setp(plt.gca(), xticklabels=[], yticklabels=[])
            if i==0 and j<5:
                plt.setp(plt.gca(), xticklabels=[])
            if i>0 and j==5:
                plt.setp(plt.gca(),  yticklabels=[])

  #  plt.figure(5)

  #  if ( not(hold)):
   #     plt.clf()
   # for j in range(0,6):
    #    print j
    #    for i in range(j+1, 6):
     #       g=get_ellipse(par.par_fid,mat1,i+3,j+3,[0,1,2] )
      #      print (j-1+(i)*5)
       #     plt.subplot(5,5,j*5+i) 
        #    plt.subplots_adjust(wspace=0.05, hspace=0.05)
         #   plt.plot(g[0], g[1])
          #  if i==j+1:
                #plt.plot(g[0],g[1])
           #     plt.xlabel(par_lab[i+3])
            #    plt.ylabel(par_lab[j+3])
           # else:
                #plt.plot(g[0],g[1])
            #    plt.setp(plt.gca(), xticklabels=[], yticklabels=[])
           
def sevenbyseven(mat, par, ind=0,hold=False):

    if ind ==0:
        fix=[1,2]
        ir=[0,3,4,5,6,7,8]
        par_vec=[par.omnuh2, par.omch2*100, par.ombh2*100, par.tau*100, par.scalar_amp,par.scalar_spectral_index*100, par.thetas*100 ]
        par_lab=['$\sum m_{\\nu}$ (eV)','$\Omega_{c}h^2$','$ \Omega_{b}h^2$','$\\tau$ ', '$A_s$', '$ n_s$', '$ 100 \Theta_s$' ]
    if ind ==1:
        fix=[0,2]
        ir=[1,3,4,5,6,7,8]
        par_vec=[par.w, par.omch2*100, par.ombh2*100, par.tau*100, par.scalar_amp,par.scalar_spectral_index*100, par.thetas*100 ]
        par_lab=['$w$','$\Omega_{c}h^2$','$ \Omega_{b}h^2$','$\\tau$ ', '$A_s$', '$ n_s$', '$100 \Theta_s$' ]
    
    if ind ==2:
        fix=[0,1]
        ir=[2,3,4,5,6,7,8]
        par_vec=[par.omk, par.omch2*100, par.ombh2*100, par.tau*100, par.scalar_amp,par.scalar_spectral_index*100, par.thetas*100 ]
        par_lab=['$\Omega_{K}$','$\Omega_{c}h^2$','$ \Omega_{b}h^2$','$\\tau$ ', '$A_s$', '$ n_s$', '$ 100 \Theta_s$'  ]
       
    mat1=mat.copy()    
    


    plt.figure(1, figsize=(10,10))

    if ( not(hold)):
        plt.clf()
    for i in range(0,7):
        print i
        for j in range(i+1,7):
            g=get_ellipse(par.par_fid,mat,ir[i],ir[j],fix )
            plt.subplot(6,6,i+1+(j-1)*6) 
            plt.subplots_adjust(wspace=0.1, hspace=0.1)
            plt.plot(g[0], g[1])
            
            deltax=nm.abs(nm.max(g[0])-nm.min(g[0]))/nm.min(g[0])
            deltay=nm.abs(nm.max(g[1])-nm.min(g[1]))/nm.min(g[1])
            print deltax
            plt.xticks(nm.arange(nm.min(g[0])-deltax*nm.min(g[0])*.01,nm.max(g[0])+deltax*nm.min(g[0])*0, deltax*nm.min(g[0])/3), rotation=90, fontsize=12)
            plt.yticks(nm.arange(nm.min(g[1])-deltay*nm.min(g[1])*.01,nm.max(g[1])+deltay*nm.min(g[1])*0, deltay*nm.min(g[1])/3), fontsize=12)
            
            if i==0 :
                plt.ylabel(par_lab[j], fontsize=15)
            if j==6:
                plt.xlabel(par_lab[i], fontsize=15)
          
            if i>0 and j<6:
                plt.setp(plt.gca(), xticklabels=[], yticklabels=[])
            if i==0 and j<6:
                plt.setp(plt.gca(), xticklabels=[])
            if i>0 and j==6:
                plt.setp(plt.gca(),  yticklabels=[])

def get_error_on_one_param(pi, F, prior):

    FiM=F.copy()
    
    
    for i in prior:
        FiM[i,i]+=1e180

    C=nm.matrix(FiM).I
    if pi == 0:
        return nm.sqrt(C[pi,pi])*94
    else:

        return nm.sqrt(C[pi,pi])



def gaussian_2d(mean, cx,cy,cxy):
    scale=1.51517
    #scale=1
    npi=2000
    theta=nm.arange(0, npi)/(npi-1.)*2*nm.pi

    x=nm.cos(theta)
    y=nm.sin(theta)

    det=cx*cy-cxy**2
    Lxx=nm.sqrt(cx)
    Lyx=cxy/Lxx
    Lyy=nm.sqrt(det)/Lxx
    

    s1=nm.sqrt(0.5* (  cx+cy- nm.sqrt( (cx-cy)**2+4*cxy**2  )   ))
    s2=nm.sqrt(0.5* (  cx+cy+ nm.sqrt( (cx-cy)**2+4*cxy**2  )   ))


    return [  Lxx*x* scale+ mean[0]      ,scale*(Lyx*x+Lyy*y  ) + mean[1]     , nm.pi*s1*s2]

def get_ellipse(mean, Fish, pi,pj, prior):

    FiM=Fish.copy()
    for i in prior:
        FiM[i,i]+=1e80

    C=nm.matrix(FiM).I
    
    if pi==0:
        mean2=[mean[pi]*93.14, mean[pj]]
        return gaussian_2d(mean2, C[pi,pi]*93.14**2, C[pj,pj], C[pi, pj]*93.14)
    elif pj ==0:
        mean2=[mean[pi], mean[pj]*93.14]
        return gaussian_2d(mean2, C[pi,pi], C[pj,pj]*93.14**2, C[pi, pj]*93.14)
    else:
        mean2=[mean[pi], mean[pj]]
        return gaussian_2d(mean2, C[pi,pi], C[pj,pj], C[pi, pj])



def create_binning_matrix2(Nbins,lint, lmax):
    """Permet de creer une matrice de binning avec deltal=1 de 2 a lint, et ensuite avec Nbins de lint a lmax"""
### CREATION DE LA MATRICE DE BINNING ###
   
    
    N_tot=Nbins+ lint-1
    B=nm.zeros((Nbins+ lint-1,lmax+1 ))
    
   
    bin1=(nm.floor(nm.linspace(lint+1, lmax+1, Nbins+1)))

    bord_droit_tab=nm.zeros((N_tot), dtype=nm.int32)
    index=0
    index2=0
    for i in range(0, N_tot+3):
        if i >= 3 and i<=lint+1:
            bord_droit_tab[index]=i
            index+=1
        if i>lint+1:
            bord_droit_tab[index]=bin1[index2+1]
            index+=1
            index2+=1

    bord_gauche_tab=nm.zeros((N_tot), dtype=nm.int32)
    index=0
    index2=0
    for i in range(0, N_tot+2):
        if i >= 2 and i<lint+1:
            bord_gauche_tab[index]=i
            index+=1
        if i>lint:
            bord_gauche_tab[index]=bin1[index2]
            index+=1
            index2+=1
    

    tab2=nm.floor((bord_droit_tab+bord_gauche_tab)/2.)
    print nm.sum(bord_droit_tab-bord_gauche_tab)
    

    for i in range(N_tot):
        B[i,bord_gauche_tab[i]:bord_droit_tab[i]]=1.0/(bord_droit_tab-bord_gauche_tab)[i]
        


    return B, bin1, tab2, N_tot


def create_binning_phi(Nbins, lmin, lint, lmax):
    """ create a linear binning    """
    return 1
    

def write_matrix( A, out_file):
    nm.savetxt(out_file, A, fmt='%1.16e')


def read_matrix(in_file):
    toto=nm.loadtxt(in_file)
    return toto




def make_noise(Noise, Beam, lmax):
    """ Create a noise array.
    Noise should be given in uk.arcmin 
    Beam should be given in arcmin
    Output is in uk."""

    n=nm.zeros(lmax+1)
    l=nm.arange(lmax+1)
    
    n=(Noise *nm.pi/180./60.)**2

    b=nm.exp(-0.5*(l*(l+1))* ( Beam*nm.pi/180./60 )**2/ nm.log(2)/8 )

    n=n/b**2
    return n


def make_multi_noise(noise_T, noise_P, theta_FWHM, lmax_cmb):


    nltt=    make_noise(0,0, lmax_cmb)
    nlpp=make_noise(0,0,lmax_cmb)


    for i in range(0,noise_T.__len__() ):
        if noise_T[i]==0:
            nltt+=0*make_noise(0, 0, lmax_cmb)
        else:
            nltt+=1.0/make_noise(noise_T[i],theta_FWHM[i] , lmax_cmb)
            nlpp+=1.0/make_noise(noise_P[i],theta_FWHM[i] , lmax_cmb)

    if noise_T[0]==0:
        nltt+=0*make_noise(0, 0, lmax_cmb)
        nlpp+=0*make_noise(0, 0, lmax_cmb)
    
    else:
        nltt=1.0/nltt
        nlpp=1.0/nlpp

    return nltt, nlpp

