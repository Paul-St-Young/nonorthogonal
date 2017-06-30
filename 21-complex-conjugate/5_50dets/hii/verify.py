#!/usr/bin/env python
import numpy as np

if __name__ == '__main__':

  hii = np.loadtxt('diag.dat')
  vmc = np.loadtxt('diff.dat')
  ndet = len(hii)
  myx = np.arange(ndet) + 1
  print vmc.shape

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('determinant index',fontsize=16)
  ax.set_ylabel('total energy relative to Hatree-Fock (ha)',fontsize=16)
  ax.plot(myx,hii,ls='-',marker='.',c='k',label='phfmol')
  ax.errorbar(myx,vmc[:,0],yerr=vmc[:,1],fmt='x',mew=1,label='vmc')
  ax.legend(loc='upper left')
  fig.tight_layout()
  fig.savefig('phfmol_vs_qmcpack.eps')
  plt.show()

# end __main__
