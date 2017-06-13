#!/usr/bin/env python
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
from convert_psir_to_psig import psig_to_psir,psir_to_psig,isosurf

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt('moR.dat')

    ikpt = 0
    imo  = 0
    moR = data[ikpt].reshape(729,10)
    # !!!! assuming cubic grid with nx points in each of the 3 dimensions
    nx = int(round(moR.shape[0]**(1./3)))

    psir = np.absolute( moR[:,imo].reshape([nx,nx,nx]) )
    print psir.min(),psir.max()

    fig = plt.figure()
    ax  = fig.add_subplot(111,projection='3d',aspect=1)
    isosurf(ax,psir,0.2)
    plt.show()
