#!/usr/bin/env python
import numpy as np

if __name__ == '__main__':

  mdvec = np.loadtxt('diff.dat')
  ref   = np.loadtxt('../1_gs16/hii/diag.dat')

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('determinant index',fontsize=16)
  ax.set_ylabel('total energy (ha)',fontsize=16)

  ax.plot(mdvec,label='hack mdvec',marker='o')
  ax.plot(ref,label='hmmt/ovmt',c='k',lw=2,ls='--')
  ax.legend(loc='upper left')
  fig.tight_layout()
  fig.savefig('single_det.eps')

  plt.show()

# end __main__
