#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
from pyscf import lib
from pyscf.pbc import gto, scf

def carbon_eigensystem(mf,detlist):

  # detlist should be checked outside of this function
  ndet = len(detlist)

  from pyscf.pbc.dft import gen_grid, numint
  coords = gen_grid.gen_uniform_grids(mf.cell)
  aoR = numint.eval_ao(mf.cell,coords)
  nao = aoR.shape[1]
  rgrid_shape = 2*np.array(mf.cell.gs)+1
  assert np.prod(rgrid_shape)==aoR.shape[0]
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
  
  if os.path.isfile(gfile) and os.path.isfile(eig_fname):
    eig_df    = pd.read_json(eig_fname).set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
    int_gvecs = np.loadtxt(gfile)
  else:
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
          'reduced_k':[0,0,0],'evalue':istate,'evector':psig}
        data.append(entry)
      # end for iorb
    # end for idet

    eig_df = pd.DataFrame(data).set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
    #eig_df.reset_index().to_json(eig_fname)
    #np.savetxt(gfile,int_gvecs)
  # end if
  return int_gvecs,eig_df
# end def carbon_eigensystem

def gen_from_detlist(detlist_fname,h5_fname):
  # get Hatree-Fock orbitals
  import sys
  sys.path.insert(0,'../1_dump_integrals/')
  from dump_integrals import run_carbon
  mf = run_carbon(verbose=3)
  detlist  = np.loadtxt(detlist_fname).view(complex)
  ndet = len(detlist)
  if ndet == 1: # single-determinant
    detlist = detlist.reshape(1,len(detlist))
  # end if
  int_gvecs,eig_df = carbon_eigensystem(mf,detlist)
  from pyscf_orbital_routines import generate_pwscf_h5
  generate_pwscf_h5(mf.cell,int_gvecs,eig_df
    ,pseudized_charge={'C':2},h5_fname=h5_fname)
  # end for
# end def gen_from_detlist

if __name__ == '__main__':

  dfname = '../2_gen_dets/out_detlist.dat'
  hfname = 'out_dets.h5'
  gen_from_detlist(dfname,hfname)

# end __main__
