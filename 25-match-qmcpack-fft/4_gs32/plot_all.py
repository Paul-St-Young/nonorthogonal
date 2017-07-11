#!/usr/bin/env python
import numpy as np

def add_phfmol_min_max(ax,phfmol_data,ewald):
  pmin = phfmol_data[:,1].min()-ewald
  pmax = phfmol_data[:,1].max()-ewald

  ax.axhline(pmin,c='k',lw=2,alpha=0.6)
  ax.axhline(pmax,c='k',lw=2,alpha=0.6)

  diff_mha = (pmax-pmin)*1000.

  ax.annotate('%d mha'%diff_mha, xy=(25,0.5*(pmin+pmax)),fontsize=24 )
  ax.arrow(15,pmax,0,(pmin-pmax+6e-3),head_width=2,head_length=0.005,fc='k',ec='k')
# end def

def add_dmc_min_max(ax,data):
  dmin = data[:,1].min()
  dmax = data[:,1].max()

  ax.axhline(dmin,c='b',lw=2,alpha=0.6)
  ax.axhline(dmax,c='b',lw=2,alpha=0.6)

  diff_mha = (dmax-dmin)*1000.

  ax.annotate('%d mha'%diff_mha, xy=(25,0.5*(dmin+dmax)-8e-3),fontsize=24 )
  ax.arrow(15,dmax,0,(dmin-dmax+6e-3),head_width=2,head_length=0.005,fc='b',ec='b')
# end def
  

if __name__ == '__main__':

  ewald = 2.695783
  phfmol_data = np.loadtxt('gen_dets/ke.dat')
  slater_data = np.loadtxt('mdets/qtot.dat')
  sj_vmc_data = np.loadtxt('dmc/sj_vmc.dat') 
  sj_dmc_data = np.loadtxt('dmc/sj_dmc.dat') 

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('number of determinants',fontsize=16)
  ax.set_ylabel('total energy (ha)',fontsize=16)

  add_phfmol_min_max(ax,phfmol_data,ewald)
  add_dmc_min_max(ax,sj_dmc_data)

  ax.plot(phfmol_data[:,0],phfmol_data[:,1]-ewald,c='k',marker='o',label='phfmol')
  ax.errorbar(slater_data[:,0],slater_data[:,1],yerr=slater_data[:,2],c='r',fmt='x--',mew=1.5,label='Slater VMC')
  ax.errorbar(sj_vmc_data[:,0],sj_vmc_data[:,1],yerr=sj_vmc_data[:,2],c='g',fmt='+--',mew=1.5,label='Slater-Jastrow VMC')
  ax.errorbar(sj_dmc_data[:,0],sj_dmc_data[:,1],yerr=sj_dmc_data[:,2],c='b',fmt='s-',label='Slater-Jastrow DMC')


  ax.legend()
  ax.set_xlim(0,51)
  fig.tight_layout()
  fig.savefig('afqmc_ci.eps')
  plt.show()
