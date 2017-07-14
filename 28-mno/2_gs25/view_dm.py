#!/usr/bin/env python

if __name__ == '__main__':
  import numpy as np
  import matplotlib.pyplot as plt

  ndorb    = 5 # 5 d orbitals
  first_mn = 12
  mn_nbas  = 34
  second_mn = first_mn + mn_nbas
  
  flist = ['1rdm_up.dat','1rdm_dn.dat','new_1rdm_up.dat','new_1rdm_dn.dat']

  fig = plt.figure()
  iplot = 0
  for fname in flist:
    iplot += 1
    mat = np.loadtxt(fname)

    ax = fig.add_subplot(2,2,iplot,aspect=1)
    ax.matshow(mat)
  # end for
    
  plt.show()
# end __main__
