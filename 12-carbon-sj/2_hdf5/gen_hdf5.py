#!/usr/bin/env python
import sys
sys.path.insert(0,'../1_eigsys/')
from carbon import build_carbon_cell

if __name__ == '__main__':

  import numpy as np
  import pandas as pd
  gvecs = np.loadtxt('../1_eigsys/gvectors.dat')
  eig_df= pd.read_json('../1_eigsys/eigensystem.json').set_index(
    ['ikpt','ispin','istate'],drop=True).sort_index()

  import h5py
  from pwscf_h5 import PwscfH5
  new = h5py.File('pwscf.pwscf.h5','w')
  ref = PwscfH5()
  cell = build_carbon_cell()
  nelecs = ref.system_from_cell(new,cell,pseudized_charge={'C':2})
  ref.create_electrons_group(new,gvecs,eig_df,nelecs)

  # transfer version info.
  new.create_dataset('application/code',data=['pyscf'])
  new.create_dataset('application/version',data=['1.4a'])
  new.create_dataset('format',data=['ES-HDF'])
  new.create_dataset('version',data=[2,1,0])

# end __main__
