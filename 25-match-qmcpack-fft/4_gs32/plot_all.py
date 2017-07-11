#!/usr/bin/env python
import numpy as np

def add_phfmol_min_max(ax,phfmol_data,ewald):
  pmin = phfmol_data[:,1].min()-ewald
  pmax = phfmol_data[:,1].max()-ewald

  ax.axhline(pmin,c='k',lw=2)
  ax.axhline(pmax,c='k',lw=2)

  diff_mha = (pmax-pmin)*1000.

  ax.annotate('%d mha'%diff_mha, xy=(25,0.5*(pmin+pmax)),fontsize=24 )
  ax.arrow(27,pmax,0,(pmin-pmax+5e-3),head_width=2,head_length=0.005,fc='k',ec='k')
# end def

def add_dmc_min_max(ax_data,ewald):
  pmin = phfmol_data[:,1].min()-ewald
  pmax = phfmol_data[:,1].max()-ewald

  ax.axhline(pmin,c='k',lw=2)
  ax.axhline(pmax,c='k',lw=2)

  diff_mha = (pmax-pmin)*1000.

  ax.annotate('%d mha'%diff_mha, xy=(25,0.5*(pmin+pmax)),fontsize=24 )
  ax.arrow(27,pmax,0,(pmin-pmax+5e-3),head_width=2,head_length=0.005,fc='k',ec='k')
# end def
  

if __name__ == '__main__':

  ewald = 2.695783
  phfmol_data = np.loadtxt('gen_dets/ke.dat')
  slater_data = np.loadtxt('mdets/qtot.dat')
  sj_vmc_data = np.loadtxt('dmc/sj_vmc.dat') 
  sj_dmc_data = np.loadtxt('dmc/sj_dmc.dat') 

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)

  add_phfmol_min_max(ax,phfmol_data,ewald)

  ax.plot(phfmol_data[:,0],phfmol_data[:,1]-ewald,c='k',marker='o',label='phfmol')
  ax.errorbar(slater_data[:,0],slater_data[:,1],yerr=slater_data[:,2],c='r',fmt='x',mew=3,label='Slater VMC')
  ax.errorbar(sj_vmc_data[:,0],sj_vmc_data[:,1],yerr=sj_vmc_data[:,2],c='g',fmt='+',mew=3,label='Slater-Jastrow VMC')
  ax.errorbar(sj_dmc_data[:,0],sj_dmc_data[:,1],yerr=sj_dmc_data[:,2],c='b',fmt='s',label='Slater-Jastrow DMC')


  ax.legend()
  ax.set_xlim(0,51)
  fig.tight_layout()
  plt.show()
