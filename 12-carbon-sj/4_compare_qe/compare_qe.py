#!/usr/bin/env python
from pwscf_h5 import PwscfH5
import numpy as np
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
from convert_psir_to_psig import isosurf

if __name__ == '__main__':

  vol = 78.7125735374 # simulation cell volume
  rgrid_shape = np.array([9,9,9])
  ikpt   = 0
  ispin  = 0
  istate = 0

  level_map = {0:0.01,1:0.015,2:0.01,3:0.015,4:0.012}
  level_map1 = {0:0.01,1:0.02,2:None,3:None,4:0.015}

  from pwscf_h5 import PwscfH5
  wf = PwscfH5()
  wf.read('../2_hdf5/pwscf.pwscf.h5')
  wf_qe = PwscfH5()
  qe_h5 = '../../10-carbon-dimer/2_orbs/pwscf.pwscf.h5'
  wf_qe.read(qe_h5)
  for istate in range(4):
    psir = wf.get_psir_from_psig(ikpt,ispin,istate,rgrid_shape)
    val = np.absolute(psir)*vol/np.prod(rgrid_shape)

    psir1 = wf_qe.get_psir_from_psig(ikpt,ispin,istate,rgrid_shape)
    val1 = np.absolute(psir1)

    overlap = np.sum(val*val1)/np.sqrt(np.sum(val*val)*np.sum(val1*val1))
    print istate,overlap

    if overlap < 0.9:
      import matplotlib.pyplot as plt
      fig = plt.figure()
      ax  = fig.add_subplot(121,projection='3d',aspect=1)
      ax1  = fig.add_subplot(122,projection='3d',aspect=1)
      ax.set_title('pyscf')
      ax1.set_title('qe')
      isosurf(ax,val,level_map[istate])
      isosurf(ax1,val1,level_map1[istate])
      plt.show()
    # end if
  # end for istate

# end __main__
