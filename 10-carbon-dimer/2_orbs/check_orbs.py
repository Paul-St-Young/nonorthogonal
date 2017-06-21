#!/usr/bin/env python
from pwscf_h5 import PwscfH5
import numpy as np
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
from convert_psir_to_psig import isosurf

if __name__ == '__main__':
  from pwscf_h5 import PwscfH5
  wf = PwscfH5()
  wf.read('pwscf.pwscf.h5')
  psir = wf.get_psir_from_psig(0,0,0,np.array([9,9,9]))
  val = np.absolute(psir)

  import matplotlib.pyplot as plt
  fig = plt.figure()
  ax  = fig.add_subplot(111,projection='3d',aspect=1)
  isosurf(ax,val)
  plt.show()
