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
import pyscf
from pyscf.pbc import gto,scf
from pyscf.pbc.scf import RHF
from pyscf import lib

def run_carbon():
  basis={
  'C':pyscf.gto.basis.parse('''
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
  alat0 = 3.6 # angstrom?
  axes  = (np.ones((3,3))-np.eye(3))*alat0/2.0
  elem  = ['C','C']
  pos   = np.array([[0,0,0],[0.25,0.25,0.25]])*alat0
  atoms = atom_text(elem,pos)
  #basis = 'gth-szv'
  #pseudo= 'gth-pade'
  fname = 'vdz.h5'

  cell = gto.M(verbose=4,a=axes,
    gs=[4,4,4],atom=atoms,basis=basis,pseudo={'C':'bfd'})#pseudo)
  cell.charge=0
  mf = RHF(cell,exxdiv=None)
  mf.max_cycle=50
  mf.chkfile=fname
  mf.conv_tol=1e-07
  mf.diis_start_cycle=1

  # run or load
  if os.path.isfile(fname):
    mf.__dict__.update(lib.chkfile.load(fname,'scf'))
  else:
    mf.scf()
  return mf
# end def run_carbon

if __name__ == '__main__':

  mf = run_carbon()
  c = mf.mo_coeff
  nmo = c.shape[1]
  h1e = reduce(np.dot, (c.T, mf.get_hcore(), c)) # 1-electron integrals
  eri = mf.with_df.ao2mo(c) # 2-electron integrals
  eri = ao2mo.restore('s8',eri,nmo) # use 8-fold symmetry of integrals of real orbitals [ij|kl] = [kj|il\ etc.
  tools.fcidump.from_integrals('fcidump.dat', 
    h1e, eri, nmo, cell.nelectron,nuc=cell.energy_nuc(), ms=0)

# end __main__
