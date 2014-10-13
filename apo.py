from __future__ import print_function

import healpy as hp
import os
import subprocess as sbp
import tempfile
import numpy as nm

def apocos(theta,sz): #sz is degree
	apo = nm.cos(theta*180/sz)/2.+.5
	return 1-nm.where(theta<sz/180.*nm.pi,apo,0)

def aposlepian(theta,sz):
  from scipy import signal
  slp = signal.slepian(1001,0.019)
  flp = nm.interp(theta,nm.linspace(0,sz/180.*nm.pi,501),slp[500:])
  return 1-nm.where(theta<sz/180.*nm.pi,flp,0)

def apogaussian(theta,fwhm):
  sg = fwhm/180*nm.pi/nm.sqrt(8*nm.log(2))
  return 1-nm.exp(-(theta/2/sg)**2)

def apogaussian_sg(theta,sg):
  ## sg = fwhm/180*nm.pi/nm.sqrt(8*nm.log(2))
  sg_rad= sg /180.*nm.pi
  return 1-nm.exp(-(theta/2/sg_rad)**2)

def compute_distance(msk,saveto=None,process_mask="./process_mask"):

  if isinstance(msk,str):
    mskfile = msk
    delmskfile=""
  else:
    mskfile = tempfile.mktemp(prefix="input_mask_",suffix=".fits")
    delmskfile = mskfile
    hp.write_map(mskfile,msk)

  delsaveto = ""
  if not saveto:
    saveto = tempfile.mktemp(prefix="output_dist_",suffix=".fits")
    delsaveto = saveto

  parfile = tempfile.mktemp(prefix="mask_",suffix=".par")
  f=open(parfile,"w")
  print("mask_file = ",mskfile,file=f)
  print("distance_file = ",saveto,file=f)
  f.close()
  sbp.call([process_mask,parfile])

  dist = hp.read_map(saveto)

  os.remove(parfile)
  if delmskfile:
    os.remove(delmskfile)
  if delsaveto:
    os.remove(delsaveto)

  return dist


  
