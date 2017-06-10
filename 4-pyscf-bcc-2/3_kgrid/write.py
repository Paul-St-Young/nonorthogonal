#!/usr/bin/env python
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
from convert_psir_to_psig import psig_to_psir,psir_to_psig,isosurf

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    import pickle

    #data = np.loadtxt('moR.dat')
    with open('data.p','r') as f:
        data = pickle.load(f)

    ikpt = 1
    abs_kpts = np.loadtxt('abs_kpts.dat')

    kvec = abs_kpts[ikpt]
    moR = data[ikpt]#.reshape(729,10)
    # !!!! assuming cubic grid with nx points in each of the 3 dimensions
    nx = int(round(moR.shape[0]**(1./3)))

    """
    alat = 3.77945227
    x,y,z= alat/nx*np.mgrid[0:nx,0:nx,0:nx]

    bloch_shift = np.exp(-1j*(kvec[0]*x+kvec[1]*y+kvec[2]*z))
    """
    psir = np.absolute( moR[:,0].reshape([nx,nx,nx]) )# * bloch_shift ) # first orbital
    print psir.min(),psir.max()

    fig = plt.figure()
    ax  = fig.add_subplot(111,projection='3d',aspect=1)
    isosurf(ax,psir,0.2)
    plt.show()
