#!/usr/bin/env python
def atom_text(elem,pos):
    assert len(elem) == len(pos)
    lines = []
    for iatom in range(len(elem)):
        mypos = pos[iatom]
        line = '%5s  %10.6f  %10.6f  %10.6f' % (elem[iatom],mypos[0],mypos[1],mypos[2])
        lines.append(line)
    atext = ';\n'.join(lines)
    return atext
# end def

import os
import numpy as np
import pyscf.gto.basis as bas
from pyscf.pbc import gto,scf
from pyscf.pbc.dft import RKS
from pyscf import lib

# !!!! hard-code basis and pseudo to avoid mistakes in later folders
def build_carbon_cell():#basis={'C':'gth-szv'},pseudo={'C':'gth-pade'}):
  #basis = {'C':'bfd-vdz'} # too many diffused functions
  basis={ # truncated basis set from Lucas
  'C':bas.parse('''
  C s
    0.205100 0.397529 
    0.409924 0.380369 
    0.819297 0.180113 
    1.637494 -0.033512 
    3.272791 -0.121499 
    6.541187 0.015176 
    13.073594 -0.000705 
   C p
    0.234064 0.302667 
    0.468003 0.289868 
    0.935757 0.210979 
    1.871016 0.112024 
    3.741035 0.054425 
    7.480076 0.021931 
   C d
    0.329486 1.000000 
   C f
    0.773485 1.000000 
   C s
  0.2 1.0
   C s
  0.6 1.0
   C p
  0.2 1.0
   C p
  0.6 1.0
  '''),
  }
  pseudo = {'C':'bfd'}
  alat0 = 3.6 # angstrom
  axes  = (np.ones((3,3))-np.eye(3))*alat0/2.0
  elem  = ['C','C']
  pos   = np.array([[0,0,0],[0.25,0.25,0.25]])*alat0
  atoms = atom_text(elem,pos)
  #basis = 'gth-szv'
  #pseudo= 'gth-pade'

  cell = gto.M(verbose=4,a=axes,
    gs=[4,4,4],atom=atoms,basis=basis,pseudo=pseudo)
  cell.charge=0
  return cell
# end def build_carbon_cell

def run_carbon(chkfile_name='vdz.h5'):
  cell = build_carbon_cell()
  mf = RKS(cell)
  mf.xc = 'pbe,pbe'
  mf.max_cycle=50
  mf.chkfile=chkfile_name
  mf.conv_tol=1e-07
  mf.diis_start_cycle=1

  # run or load
  if os.path.isfile(chkfile_name):
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  return mf
# end def run_carbon

def ao_on_grid(cell):
  from pyscf.pbc.dft import gen_grid,numint
  coords = gen_grid.gen_uniform_grids(cell)
  aoR    = numint.eval_ao(cell,coords)
  return aoR
# end def ao_on_grid

def mo_coeff_to_psig(mo_coeff,aoR,cell_gs,cell_vol,int_gvecs=None):
  """
   Inputs:
     mo_coeff: molecular orbital in AO basis, each column is an MO, shape (nao,nmo)
     aoR: atomic orbitals on a real-space grid, each column is an AO, shape (ngrid,nao)
     cell_gs: 2*cell_gs+1 should be the shape of real-space grid (e.g. (5,5,5))
     cell_vol: cell volume, used for FFT normalization
     int_gvecs: specify the order of plane-waves using reciprocal lattice points
   Outputs:
       3. plane-wave coefficients representing the MOs, shape (ngrid,nmo)
  """
  # provide the order of reciprocal lattice vectors to skip
  if int_gvecs is None: # use internal order
    nx,ny,nz = cell_gs
    from itertools import product
    int_gvecs = np.array([gvec for gvec in product(
      range(-nx,nx+1),range(-ny,ny+1),range(-nz,nz+1))],dtype=int)
  else:
    assert (int_gvecs.dtype is int)
  # end if
  npw = len(int_gvecs) # number of plane waves 

  # put molecular orbitals on real-space grid
  moR = np.dot(aoR,mo_coeff) 
  nao,nmo = moR.shape
  rgrid_shape = 2*np.array(cell_gs)+1
  assert nao == np.prod(rgrid_shape)

  # sort MOs from lowest to highest energy
  # for each MO, FFT to get psig
  psig = np.zeros([nmo,npw,2]) # store real & complex
  for istate in range(nmo):
    # fill real-space FFT grid
    rgrid = moR[:,istate].reshape(rgrid_shape)
    # get plane-wave coefficients (on reciprocal-space FFT grid)
    moG   = np.fft.fftn(rgrid)/np.prod(rgrid_shape)*cell_vol
    # transfer plane-wave coefficients to psig in specified order
    for igvec in range(npw):
      comp_val = moG[tuple(int_gvecs[igvec])]
      psig[istate,igvec,:] = comp_val.real,comp_val.imag
    # end for igvec
  # end for istate
  return int_gvecs,psig
# end def mo_coeff_to_psig

if __name__ == '__main__':

  # run DFT
  mf = run_carbon()

  # put AO on real-space grid
  aoR = ao_on_grid(mf.cell)
  
  # put MOs in file for PwscfH5 to read
  import pandas as pd
  gfile = 'gvectors.dat'
  fname = 'eigensystem.json'
  if os.path.isfile(fname):
    gvecs = np.loadtxt(gfile)
    eig_df = pd.read_json(fname).set_index(
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
    eig_df.to_json(fname)
    np.savetxt(gfile,gvecs)
  # end if
  """
  c = mf.mo_coeff
  nmo = c.shape[1]
  h1e = reduce(np.dot, (c.T, mf.get_hcore(), c)) # 1-electron integrals
  eri = mf.with_df.ao2mo(c) # 2-electron integrals
  eri = ao2mo.restore('s8',eri,nmo) # use 8-fold symmetry of integrals of real orbitals [ij|kl] = [kj|il\ etc.
  tools.fcidump.from_integrals('fcidump.dat', 
    h1e, eri, nmo, cell.nelectron,nuc=cell.energy_nuc(), ms=0)
  """

# end __main__
