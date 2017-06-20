#!/usr/bin/env python
import os
import numpy as np

def build_bcc2_cell(verbose=4):
  #import pyscf.pbc.gto as pbcgto
  alat  = 3.77945227
  basis = 'cc-pVDZ'
  #ke    = 20.

  # define system
  axes = alat*np.eye(3)
  atom_text = 'H 0 0 0; H %f %f %f' % tuple(alat*np.array([0.5,0.5,0.5]))
  #cell = pbcgto.Cell()
  #cell.build(
  #  a         = axes,
  #  atom      = atom_text,
  #  unit      = 'B',
  #  basis     = basis,
  #  ke_cutoff = ke,
  #  verbose   = verbose
  #)
  from pyscf.pbc import gto
  cell = gto.M(verbose=verbose,a=axes,gs=[4,4,4]
    ,atom=atom_text,basis=basis,unit='bohr')

  return cell
# end def build_bcc2_cell

def run_bcc2(chkfile_name='vdz.h5',verbose=4):
  from pyscf.pbc.scf import RHF
  cell = build_bcc2_cell(verbose=verbose)
  mf   = RHF(cell,exxdiv=None)
  mf.max_cycle = 50
  mf.conv_tol  = 1e-7
  mf.diis_start_cycle = 1
  mf.chkfile   = chkfile_name
  # run or load
  if os.path.isfile(chkfile_name):
    from pyscf import lib
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  return mf
# end def run_bcc2

import sys
sys.path.insert(0,'../../12-carbon-sj/1_eigsys')
from carbon import ao_on_grid,mo_coeff_to_psig
def save_eigensystem(mf,gfile='gvectors.dat',eig_fname='eigensystem.json'):
  # put AO on real-space grid
  aoR = ao_on_grid(mf.cell)

  # put MOs in file for PwscfH5 to read
  import pandas as pd
  if os.path.isfile(eig_fname):
    gvecs = np.loadtxt(gfile)
    eig_df = pd.read_json(eig_fname).set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
  else:
    data = []
    ikpt  = 0 # gamma-point calculation
    ispin = 0 # restricted (same orbitals for up and down electrons)
    # get MOs in plane-wave basis
    gvecs,psig = mo_coeff_to_psig(mf.mo_coeff,aoR,mf.cell.gs,mf.cell.vol)
    nstate,npw,ncomp = psig.shape
    for istate in range(nstate):
      entry = {'ikpt':ikpt,'ispin':ispin,'istate':istate,
        'reduced_k':mf.kpt,'evalue':mf.mo_energy[istate],'evector':psig[istate,:,:]}
      data.append(entry)
    # end for istate
    eig_df = pd.DataFrame(data)
    eig_df.to_json(eig_fname)
    np.savetxt(gfile,gvecs)
  # end if
# end def

if __name__ == '__main__':
  mean_field_object = run_bcc2()
  save_eigensystem(mean_field_object)
# end __main__
