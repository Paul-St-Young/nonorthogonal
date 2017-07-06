#!/usr/bin/env python
import numpy as np

if __name__ == '__main__':

  mdvec = np.loadtxt('diff.dat')
  ref   = np.loadtxt('../1_gs16/hii/diag.dat')
  vmc   = np.loadtxt('../1_gs16/hii/diff.dat')
  ndet = len(mdvec)

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('determinant index',fontsize=16)
  ax.set_ylabel('total energy (ha)',fontsize=16)
  ax.set_xlim(-0.5,ndet-0.5)

  myx = range(ndet)
  ax.plot(myx,mdvec,label='hack mdvec',marker='o')
  ax.plot(myx,ref,label='hmmt/ovmt',c='k',lw=2,ls='--')
  ax.errorbar(myx,vmc[:,0],yerr=vmc[:,1],fmt='x',label='vmc',ms=10,mew=1)
  ax.legend(loc='upper left')
  fig.tight_layout()
  #fig.savefig('single_det.eps')

  plt.show()

# end __main__
