#!/usr/bin/env python
import numpy as np
import pandas as pd
import sys
sys.path.insert(0,'../2_system_from_cell')
from system_from_cell import bcc2
sys.path.insert(0,'../3_gvec')
from ecut import gvectors_within_cutoff

if __name__ == '__main__':

    # build simulation cell
    alat  = 3.77945227
    ke    = 20.
    cell  = bcc2(alat,ke=ke)

    # get axes
    axes = cell.a
    recip_axes = 2*np.pi*np.linalg.inv(axes)
    # get gvectors
    int_gvecs  = gvectors_within_cutoff(cell.gs,ke,recip_axes)
    #sel_neg = (int_gvecs<0)
    #nx,ny,nz = 2*cell.gs+1
    #assert (nx==ny)&(ny==nz)
    #int_gvecs[sel_neg] += nx
    
    # write pwscf.h5
    import h5py
    from pwscf_h5 import PwscfH5
    new = h5py.File('pwscf.pwscf.h5','w')
    ref = PwscfH5()
    ref.system_from_cell(new,cell)

    eig_df = pd.read_json('../4_eigensystem/mydf.json').set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
    ref.create_electrons_group(new,int_gvecs,eig_df)
    # transfer orbital info.
    new.create_dataset('electrons/number_of_electrons',data=[1,1])
    new.create_dataset('electrons/number_of_kpoints',data=[8])
    new.create_dataset('electrons/number_of_spins',data=[1])
    new.create_dataset('electrons/psi_r_is_complex',data=[1])

    # transfer version info.
    new.create_dataset('application/code',data=['pyscf'])
    new.create_dataset('application/version',data=['69d4b826c01950437f8e16663120942d5709c5e3'])
    new.create_dataset('format',data=['ES-HDF'])
    new.create_dataset('version',data=[2,1,0])

# end __main__
