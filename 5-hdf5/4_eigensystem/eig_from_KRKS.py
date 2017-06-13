#!/usr/bin/env python
import os
import numpy as np
from pyscf.pbc import gto as pbcgto
from pyscf.pbc import dft as pbcdft

import sys
sys.path.insert(0,'../2_system_from_cell')
from system_from_cell import bcc2
sys.path.insert(0,'../3_gvec')
from ecut import gvectors_within_cutoff

def eigensystem_from_pyscf(kmf,rkpts,abs_kpts,axes):
    # axes only used to get volumn
    vol = np.dot(np.cross(axes[0],axes[1]),axes[2])
    nkpt  = len(rkpts)
    nspin = 1 # !!!! assume same orbitals for up and down e-
    data  = []
    for ikpt in range(nkpt):
      for ispin in range(nspin): 
        evalues,evectors = kmf.get_bands(abs_kpts[ikpt])
        nstate = len(evalues)
        # normalize eigenvectors
        nevecs = evectors/np.linalg.norm(evectors,axis=0)
        moR = np.dot(aoR,nevecs).real # put molecular orbitals on real-space grid
        for istate in range(nstate):
          rgrid = moR[:,istate].reshape(rgrid_shape)
          # get plane-wave coefficients
          moG   = np.fft.fftn(rgrid)/np.prod(rgrid_shape)*vol  
          psig  = np.zeros([len(int_gvecs),2]) # store real & complex
          for igvec in range(len(int_gvecs)):
            comp_val = moG[tuple(int_gvecs[igvec])]
            psig[igvec,:] = comp_val.real,comp_val.imag
          # end for igvec
          entry = {'ikpt':ikpt,'ispin':ispin,'istate':istate,
            'reduced_k':rkpts[ikpt],'evalue':evalues[istate],'evector':psig}
          data.append(entry)
        # end for istate
      # end for ispin
    # end for ikpt
    import pandas as pd
    df = pd.DataFrame(data).set_index(['ikpt','ispin','istate'],drop=True)
    return df
# end def eigensystem_from_pyscf

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

    # put atomic orbitals on real-space grid
    from pyscf.pbc.dft import gen_grid,numint
    coords = gen_grid.gen_uniform_grids(cell)
    aoR = numint.eval_ao(cell,coords)
    rgrid_shape = cell.gs*2+1 # shape of real-space grid
    assert np.prod(rgrid_shape) == aoR.shape[0]

    # build simulation
    chkfile  = 'pdvz.chk'
    abs_kpts = cell.make_kpts([2,2,2])        # kpoints in 2pi/alat units
    kmf = pbcdft.KRKS(cell,abs_kpts)
    if os.path.isfile(chkfile):
        from pyscf.pbc import lib
        kmf.__dict__.update(lib.chkfile.load(chkfile,'scf'))
    else:
        kmf.xc = 'lda,lda'
        kmf.verbose = 3
        kmf.chkfile = chkfile
        kmf.scf()
    # end if

    rkpts    = cell.get_scaled_kpts(abs_kpts) # kpoints in crystal units
    df = eigensystem_from_pyscf(kmf,rkpts,abs_kpts,axes)
    df.reset_index().to_json('mydf.json')

    """
    import h5py
    from pwscf_h5 import PwscfH5
    new = h5py.File('pwscf.pwscf.h5','w')
    ref = PwscfH5()
    ref.system_from_cell(new,cell)
    """

# end __main__
