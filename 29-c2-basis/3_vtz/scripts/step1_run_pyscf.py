import numpy as np

def build_carbon_cell(gs,verbose=4):
  from pyscf.pbc import gto
  import pyscf.gto.basis as bas
  basis={ # truncated BFD triple-zeta basis set
  'C':bas.parse('''
C s
 0.205100     0.397529
 0.409924     0.380369
 0.819297     0.180113
 1.637494     -0.033512
 3.272791     -0.121499
 6.541187     0.015176
 13.073594     -0.000705
C s
 0.921552     1.000000
C s
 0.132800     1.000000
C p
 0.234064     0.302667
 0.468003     0.289868
 0.935757     0.210979
 1.871016     0.112024
 3.741035     0.054425
 7.480076     0.021931
C p
 0.126772     1.000000
C p
 0.376742     1.000000
C d
 0.329486     1.000000
C d
 1.141611     1.000000
C f
 0.773485     1.000000
  '''),
  }
  pseudo = {'C':'bfd'}
  alat0  = 3.6 # angstrom

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

def run_carbon(gs,verbose=4,chkfile_name='vdz.h5'):
  import os
  from pyscf.pbc.scf import RHF
  cell = build_carbon_cell(gs,verbose=verbose)
  mf = RHF(cell,exxdiv=None)
  mf.max_cycle=50
  mf.chkfile=chkfile_name
  mf.conv_tol=1e-7
  mf.diis_start_cycle=1
  #mf.direct_scf_tol=1e-7

  # run or load
  if os.path.isfile(chkfile_name):
    from pyscf import lib
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  return mf
# end def run_carbon

def step1_run_pyscf(gs):
  # run pyscf and extract Kohn-Sham eigensystem
  # ================================================
  mf = run_carbon(gs,verbose=3)
  return mf
