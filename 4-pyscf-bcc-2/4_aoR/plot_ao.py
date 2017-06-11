#!/usr/bin/env python
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
from convert_psir_to_psig import psig_to_psir,psir_to_psig,isosurf

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    aoR = np.loadtxt('aoR.dat')
    nx  = int(round(aoR.shape[0]**(1./3)))

    fig = plt.figure()

    for iao in range(aoR.shape[1]):
      psir = np.absolute( aoR[:,iao].reshape([nx,nx,nx]) )
      iax = (iao+1)%10
      if iao < 5: # shift grid for basis centered on the corner atom
          psir = np.fft.fftshift(psir)
      ax  = fig.add_subplot('25%d'%iax,projection='3d',aspect=1)
      isosurf(ax,psir)
    # end for
    plt.show()

# end __main__
