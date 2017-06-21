#!/usr/bin/env python
import os

if __name__ == '__main__':
  pyscf_h5_fname = '../2_hdf5/pwscf.pwscf.h5'
  qe_h5_fname    = '../../10-carbon-dimer/1_nexus/gamma/ecut20/nscf/pwscf_output/pwscf.pwscf.h5'

  # open reference files for reading
  from pwscf_h5 import PwscfH5
  wf = PwscfH5()
  wf.read(pyscf_h5_fname)

  wf_qe = PwscfH5()
  wf_qe.read(qe_h5_fname)

  # start new hdf5 file and transfer data
  import h5py
  new = h5py.File('test.h5','w')
  new.create_group('electrons')
  for name in wf.fp.keys():
    if name != 'electrons':
      # copy from pyscf
      wf.fp.copy(name,new)
    else: # copy from qe
      for name1 in wf_qe.fp[name].keys():
        if name1 == 'density':
          continue # skip density section
        wf_qe.fp.copy(os.path.join(name,name1),new[name])
      # end for
    # end if
  # end for
  new.close()
# end __main__

