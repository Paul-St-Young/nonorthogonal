#!/usr/bin/env python
import os
import numpy as np
import pandas as pd

import sys
sys.path.insert(0,'../../6-verify_hdf5/1_psig_to_moR')
from psig_to_moR import simulate_bcc2

def no_ecut_from_pyscf(kmf,cell,kgrid):

    # generate gvectors (no ke_cutoff)
    axes = cell.a
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
        istate = ikpt # !!!! fake all k-points into gamma as different states

        # build e^ikr at this kpoint
        kvec = abs_kpts[ikpt]
        eik_dot_r = np.fft.fftshift( np.zeros(rgrid_shape,dtype=complex) )
        for ix,iy,iz in product(range(-nx,nx+1),range(-ny,ny+1),range(-nz,nz+1)):
          rvec = np.dot([ix,iy,iz],axes)
          eik_dot_r[ix,iy,iz] = np.exp(1j*np.dot(kvec,rvec))

        # get eigensystem
        evectors = kmf.mo_coeff[ikpt]#.get_bands(abs_kpts[ikpt])
        assert evectors.shape== (10,10)
        moR = np.dot(aoR,evectors).real # put molecular orbitals on real-space grid

        rgrid = moR[:,0].reshape(rgrid_shape) # !!!! use the occupied state at each kpoint
        # get plane-wave coefficients
        moG   = np.fft.fftn(rgrid)/np.prod(rgrid_shape)*cell.vol * eik_dot_r
        psig  = np.zeros([len(int_gvecs),2]) # store real & complex
        for igvec in range(len(int_gvecs)):
          comp_val = moG[tuple(int_gvecs[igvec])]
          psig[igvec,:] = comp_val.real,comp_val.imag
        # end for igvec
        entry = {'ikpt':0,'ispin':ispin,'istate':istate, # !!!! store at gamma
          'reduced_k':rkpts[0],'evalue':istate,'evector':psig}
        data.append(entry)
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
  nelecs = ref.system_from_cell(new,cell)
  ref.create_electrons_group(new,int_gvecs,eig_df,nelecs)

  # transfer version info.
  new.create_dataset('application/code',data=['pyscf'])
  new.create_dataset('application/version',data=['69d4b826c01950437f8e16663120942d5709c5e3'])
  new.create_dataset('format',data=['ES-HDF'])
  new.create_dataset('version',data=[2,1,0])
# end def main

if __name__ == '__main__':
  main()
