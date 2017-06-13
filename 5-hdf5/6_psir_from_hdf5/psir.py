#!/usr/bin/env python
import numpy as np
from pwscf_h5 import PwscfH5
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert/')
from convert_psir_to_psig import isosurf,psig_to_psir

if __name__ == '__main__':

    fhandle = PwscfH5()
    fhandle.read('../5_pyscf_hdf5/pwscf.pwscf.h5')

    loc = {'ikpt':0,'ispin':0,'istate':0}
    psig_arr = fhandle.psig(**loc)
    psig = psig_arr[:,0] + 1j*psig_arr[:,1]

    gvec = fhandle.get('gvectors')
    psir = psig_to_psir(psig,gvec,(9,9,9))
    val  = np.absolute(psir)
    print val.min(),val.max()

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax  = fig.add_subplot(111,projection='3d',aspect=1)
    isosurf(ax,val)
    plt.show()
