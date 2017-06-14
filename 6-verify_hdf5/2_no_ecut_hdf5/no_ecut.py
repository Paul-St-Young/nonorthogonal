#!/usr/bin/env python
import os
import numpy as np
import pandas as pd

import sys
sys.path.insert(0,'../1_psig_to_moR')
from psig_to_moR import simulate_bcc2

def no_ecut_from_pyscf(kmf,cell,kgrid):

    # generate gvectors (no ke_cutoff)
    nx,ny,nz = cell.gs
    from itertools import product
    int_gvecs = np.array([gvec for gvec in product(
      range(-nx,nx+1),range(-ny,ny+1),range(-nz,nz+1))])

    # put AO on real-space grid
    from pyscf.pbc.dft import gen_grid,numint
    coords = gen_grid.gen_uniform_grids(cell)
    aoR    = numint.eval_ao(cell,coords)
    rgrid_shape = 2*cell.gs+1
    assert np.prod(rgrid_shape)==aoR.shape[0]

    # go through desired kpoints
    abs_kpts = cell.make_kpts(kgrid)
    rkpts    = cell.get_scaled_kpts(abs_kpts)
    nkpt  = len(rkpts)
    nspin = 1 # !!!! assume same orbitals for up and down e-
    data  = []
    for ikpt in range(nkpt):
      for ispin in range(nspin):
        # get eigensystem
        evalues  = kmf.mo_energy[ikpt]
        evectors = kmf.mo_coeff[ikpt]#.get_bands(abs_kpts[ikpt])
        assert evalues.shape == (10,)
        assert evectors.shape== (10,10)
        nstate = len(evalues)
        moR = np.dot(aoR,evectors).real # put molecular orbitals on real-space grid
        for istate in range(nstate):
          rgrid = moR[:,istate].reshape(rgrid_shape)
          # get plane-wave coefficients
          moG   = np.fft.fftn(rgrid)/np.prod(rgrid_shape)*cell.vol
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
    return df,int_gvecs
# end def no_ecut_from_pyscf

def main():
  # run DFT
  kgrid = [2,2,2]
  cell,kmf = simulate_bcc2(
    alat=3.77945227, basis='cc-pVDZ',chkfile='pvdz.h5',ke=20.,kgrid=kgrid)

  # read eigensystem
  fname = 'eigensystem.json'
  gfile = 'gvectors.dat'
  if os.path.isfile(fname):
    eig_df = pd.read_json(fname).set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
    int_gvecs = np.loadtxt(gfile)
  else:
    eig_df,int_gvecs = no_ecut_from_pyscf(kmf,cell,kgrid)
    eig_df.reset_index().to_json(fname)
    np.savetxt(gfile,int_gvecs)
  # end if

  # write pwscf.h5
  import h5py
  from pwscf_h5 import PwscfH5
  new = h5py.File('pwscf.pwscf.h5','w')
  ref = PwscfH5()
  ref.system_from_cell(new,cell)

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
# end def main

if __name__ == '__main__':
  main()
