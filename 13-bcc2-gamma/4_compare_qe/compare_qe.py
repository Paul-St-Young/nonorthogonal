#!/usr/bin/env python
from pwscf_h5 import PwscfH5
import numpy as np
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
from convert_psir_to_psig import isosurf

if __name__ == '__main__':

  vol = 53.9866768459 # simulation cell volume in bohr^3/vol

  rgrid_shape = np.array([9,9,9])
  ikpt   = 0
  ispin  = 0
  istate = 0

  from pwscf_h5 import PwscfH5
  wf = PwscfH5()
  wf.read('../2_hdf5/pwscf.pwscf.h5')
  wf_qe = PwscfH5()
  qe_h5 = '../../3-psir2hdf5/small_convert/pwscf.pwscf.h5'
  wf_qe.read(qe_h5)
  for istate in range(1):
    psir = wf.get_psir_from_psig(ikpt,ispin,istate,rgrid_shape)
    val = np.absolute(psir)*vol/np.prod(rgrid_shape)

    psir1 = wf_qe.get_psir_from_psig(ikpt,ispin,istate,rgrid_shape)
    val1 = np.absolute(psir1)

    overlap = np.sum(val*val1)/np.sqrt(np.sum(val*val)*np.sum(val1*val1))
    print istate,overlap

    if overlap < 0.7:
      import matplotlib.pyplot as plt
      fig = plt.figure()
      ax  = fig.add_subplot(121,projection='3d',aspect=1)
      ax1  = fig.add_subplot(122,projection='3d',aspect=1)
      ax.set_title('pyscf')
      ax1.set_title('qe')
      isosurf(ax,val,0.011)
      isosurf(ax1,val1,0.02)
      plt.show()
    # end if
  # end for istate

# end __main__
