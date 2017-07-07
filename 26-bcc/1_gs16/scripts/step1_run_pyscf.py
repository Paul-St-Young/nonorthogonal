import numpy as np

def build_bcc2_cell(gs,verbose=4):
  from pyscf.pbc import gto
  basis  = {'H':'gth-dzv'}

  axes  = np.array(
      [[ 8.,  0.,  0.],
       [ 0.,  8.,  0.],
       [ 0.,  0.,  8.]])
  pos   = np.array(
      [[ 0.,  0.,  0.],
       [ 2.,  2.,  2.],
       [ 4.,  0.,  0.],
       [ 6.,  2.,  2.],
       [ 0.,  4.,  0.],
       [ 2.,  6.,  2.],
       [ 4.,  4.,  0.],
       [ 6.,  6.,  2.],
       [ 0.,  0.,  4.],
       [ 2.,  2.,  6.],
       [ 4.,  0.,  4.],
       [ 6.,  2.,  6.],
       [ 0.,  4.,  4.],
       [ 2.,  6.,  6.],
       [ 4.,  4.,  4.],
       [ 6.,  6.,  6.]])
  elem = ['H'] * len(pos)
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
  atomt = atom_text(elem,pos)
  cell = gto.M(verbose=verbose,a=axes,unit='B',
    gs=gs,atom=atomt,basis=basis)
  return cell
# end def 

def run_bcc2(gs,verbose=4,chkfile_name='dzv.h5'):
  import os
  from pyscf.pbc.scf import RHF
  cell = build_bcc2_cell(gs,verbose=verbose)
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

def step1_run_pyscf(gs):
  # run pyscf and extract Kohn-Sham eigensystem
  # ================================================
  mf = run_bcc2(gs,verbose=3)
  return mf
