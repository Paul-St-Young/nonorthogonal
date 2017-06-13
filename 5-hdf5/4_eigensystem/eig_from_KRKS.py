#!/usr/bin/env python
import numpy as np
from pyscf.pbc import gto as pbcgto
from pyscf.pbc import dft as pbcdft

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

    # build simulation
    #abs_kpts = cell.make_kpts([2,2,2])        # kpoints in 2pi/alat units
    #kmf = pbcdft.KRKS(cell,abs_kpts)
    #kmf.xc = 'lda,lda'
    #kmf.verbose = 3
    #kmf.scf()

    #rkpts    = cell.get_scaled_kpts(abs_kpts) # kpoints in crystal units

    # get axes
    axes = cell.a
    recip_axes = 2*np.pi*np.linalg.inv(axes)
    # get gvectors
    int_gvecs  = gvectors_within_cutoff(cell.gs,ke,recip_axes)

    from pyscf.pbc.dft import gen_grid,numint
    coords = gen_grid.gen_uniform_grids(cell)
    aoR = numint.eval_ao(cell,coords)
    print aoR.shape,int_gvecs.shape

    """
    import h5py
    from pwscf_h5 import PwscfH5
    new = h5py.File('pwscf.pwscf.h5','w')
    ref = PwscfH5()
    ref.system_from_cell(new,cell)
    """

# end __main__
