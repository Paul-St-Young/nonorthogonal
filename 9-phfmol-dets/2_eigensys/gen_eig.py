#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
from pyscf import lib
from pyscf.pbc import gto, scf

import sys
sys.path.insert(0,'../1_ref/first/')
from carbon import run_carbon

def main():
    mf = run_carbon()

    from pyscf.pbc.dft import gen_grid, numint
    coords = gen_grid.gen_uniform_grids(mf.cell)
    aoR = numint.eval_ao(mf.cell,coords)
    nao = aoR.shape[1]
    rgrid_shape = 2*np.array(mf.cell.gs)+1
    assert np.prod(rgrid_shape)==aoR.shape[0]

    ci_coeff = np.loadtxt('../1_ref/ci_coeff.dat').view(complex)
    detlist  = np.loadtxt('../1_ref/detlist.dat').view(complex)
    ndet = len(ci_coeff)
    if ndet == 1: # single-determinant
      detlist = detlist.reshape(1,len(detlist))
    # end if
    assert detlist.shape[0] == ndet
    assert detlist.shape[1] == nao*nao

    eig_fname = 'eigsys.json'
    gfile     = 'gvectors.dat'

    # generate gvectors (no ke_cutoff)
    nx,ny,nz = mf.cell.gs
    from itertools import product
    int_gvecs = np.array([gvec for gvec in product(
      range(-nx,nx+1),range(-ny,ny+1),range(-nz,nz+1))])
    npw = len(int_gvecs)

    # turn detlist into a dataframe containing the eigensystem
    nfill = 4 # !!!! hard-code 4 filled orbitals
    ikpt=ispin=0 # only do this for RHF wavefunction at Gamma
    
    data = []
    for idet in range(ndet):
      det = detlist[idet].reshape(nao,nao)
      mo_coeff = np.dot(mf.mo_coeff,det)
      moR = np.dot(aoR,mo_coeff)
      for iorb in range(nfill):
        rgrid = moR[:,iorb].reshape(rgrid_shape)
        moG   = np.fft.fftn(rgrid)/np.prod(rgrid_shape)*mf.cell.vol
        psig  = np.zeros([npw,2]) # store real & complex
        for igvec in range(npw):
          comp_val = moG[tuple(int_gvecs[igvec])]
          psig[igvec,:] = comp_val.real,comp_val.imag
        # end for igvec
        istate = idet*nfill+iorb
        entry  = {'ikpt':ikpt,'ispin':ispin,'istate':istate,
          'reduced_k':[0,0,0],'evalue':iorb,'evector':psig}
        data.append(entry)
      # end for iorb
    # end for idet

    df = pd.DataFrame(data)
    df.to_json(eig_fname)
    np.savetxt(gfile,int_gvecs)

    eig_df    = pd.read_json(eig_fname).set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
    int_gvecs = np.loadtxt(gfile)
    # write pwscf.h5
    import h5py
    from pwscf_h5 import PwscfH5
    new = h5py.File('pwscf.pwscf.h5','w')
    ref = PwscfH5()
    nelecs = ref.system_from_cell(new,mf.cell,pseudized_charge={'C':2})
    ref.create_electrons_group(new,int_gvecs,eig_df,nelecs)

    # transfer version info.
    new.create_dataset('application/code',data=['pyscf'])
    new.create_dataset('application/version',data=['69d4b826c01950437f8e16663120942d5709c5e3'])
    new.create_dataset('format',data=['ES-HDF'])
    new.create_dataset('version',data=[2,1,0])
# end def main

if __name__ == '__main__':

    main()

# end __main__
