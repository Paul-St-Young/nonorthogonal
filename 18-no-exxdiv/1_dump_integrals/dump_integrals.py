#!/usr/bin/env python
import numpy as np

def build_carbon_cell(verbose=4):
  from pyscf.pbc import gto
  import pyscf.gto.basis as bas
  basis={ # truncated BFD double-zeta basis set from Lucas K. Wagner
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
  alat0  = 3.6 # angstrom
  gs     = [16,16,16]

  # convert elem,pos to text representation
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
  axes  = (np.ones((3,3))-np.eye(3))*alat0/2.0
  elem  = ['C','C']
  pos   = np.array([[0,0,0],[0.25,0.25,0.25]])*alat0
  atoms = atom_text(elem,pos)
  cell = gto.M(verbose=verbose,a=axes,
    gs=gs,atom=atoms,basis=basis,pseudo=pseudo)
  return cell
# end def 

def run_carbon(verbose=4,chkfile_name='vdz.h5'):
  import os
  from pyscf.pbc.scf import RHF
  cell = build_carbon_cell(verbose=verbose)
  mf = RHF(cell,exxdiv=None)
  mf.max_cycle=50
  mf.chkfile=chkfile_name
  mf.conv_tol=1e-07
  mf.diis_start_cycle=1

  # run or load
  if os.path.isfile(chkfile_name):
    from pyscf import lib
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  return mf
# end def run_carbon

if __name__ == '__main__':

  # run pyscf and extract Kohn-Sham eigensystem
  # ================================================
  print(datetime.now())
  print('performing Hatree-Fock calculation')
  mf = run_carbon(verbose=3)
  print(mf.e_tot)
  print(datetime.now())

  import os
  if not os.path.isfile('fcidump.dat'):
    from datetime import datetime
    from pyscf import tools, ao2mo
    from functools import reduce
    cell = mf.cell
    c = mf.mo_coeff
    nmo = c.shape[1]
    print(datetime.now())
    print('evaluating 1-electron integrals')
    h1e = reduce(np.dot, (c.T, mf.get_hcore(), c)) # 1-electron integrals
    print(datetime.now())
    print('evaluating 2-electron integrals')
    eri = mf.with_df.ao2mo(c) # 2-electron integrals
    print(datetime.now())
    print('restore 8-fold symmetry of the integral table')
    eri = ao2mo.restore('s8',eri,nmo) # use 8-fold symmetry of integrals of real orbitals [ij|kl] = [kj|il\ etc.
    print(datetime.now())
    print('dumping the integral table')
    tools.fcidump.from_integrals('fcidump.dat',
      h1e, eri, nmo, cell.nelectron,nuc=cell.energy_nuc(), ms=0)
    print(datetime.now())
  # end if

# end __main__
