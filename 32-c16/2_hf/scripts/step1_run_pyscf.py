import numpy as np

def bfd_basis(nz=2):
  """ return a truncated BFD basis set for carbon
  nz = 2 double-zeta 
     = 3 triple-zeta etc. 
  """
  from pyscf.gto.basis import parse
  if nz == 2:
    basis={ # truncated BFD double-zeta basis set
    'C':parse('''
  C s 
   0.205100     0.397529
   0.409924     0.380369
   0.819297     0.180113
   1.637494     -0.033512
   3.272791     -0.121499
   6.541187     0.015176
   13.073594     -0.000705
  C s
   0.127852     1.000000
  C p
   0.234064     0.302667
   0.468003     0.289868
   0.935757     0.210979
   1.871016     0.112024
   3.741035     0.054425
   7.480076     0.021931
  C d
   0.561160     1.000000
    '''),
    }
  elif nz == 3:
    basis={ # truncated BFD triple-zeta basis set
    'C':parse('''
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
  elif nz == 4:
    basis={ # truncated BFD quadruple-zeta basis set
    'C':parse('''
  C s
   0.205100     0.397529
   0.409924     0.380369
   0.819297     0.180113
   1.637494     -0.033512
   3.272791     -0.121499
   6.541187     0.015176
   13.073594     -0.000705
  C s
   0.846879     1.000000
  C s
   0.269659     1.000000
  C p
   0.234064     0.302667
   0.468003     0.289868
   0.935757     0.210979
   1.871016     0.112024
   3.741035     0.054425
   7.480076     0.021931
  C p
   0.313254     1.000000
  C p
   0.804681     1.000000
  C d
   0.240171     1.000000
  C d
   0.684884     1.000000
  C d
   2.013760     1.000000
  C f
   0.457302     1.000000
  C f
   1.324930     1.000000
  C g
   1.034180     1.000000
    '''),
    }
  else:
    raise RuntimeError('nz must be one of 2,3,4')
  # end if
  return basis
# end def bfd_basis

def build_cell(xsf_file,gs,nz=2,verbose=4):
  from pyscf.pbc import gto
  from nexus import Structure # use nexus.Structure to read xsf
  from pyscf_orbital_routines import atom_text # take ['C','C'],[[0,0,0],[0.5,0.5,0.5]] to make 'C 0 0 0\n C 0.5 0.5 0.5'

  # define Hamiltonian and atomic basis set
  pseudo = {'C':'bfd'}
  basis  = bfd_basis(nz=nz)

  # read structure from xsf file
  struct = Structure()
  struct.read(xsf_file)
  struct.change_units('A') # use angstrom units for pyscf
  axes  = struct.axes
  pos   = struct.pos
  elem  = struct.elem

  # construct simulation cell
  cell = gto.M(
    verbose = verbose,
    gs      = gs,
    a       = axes,
    atom    = atom_text(elem,pos),
    basis   = basis,
    pseudo  = pseudo,
    unit    = 'angstrom'
  )
  return cell
# end def build_cell

def run_pyscf(xsf_file,gs,nz=2,verbose=4,chkfile_name='bfd.h5',exxdiv=None):
  import os
  from pyscf.pbc import scf

  cell = build_cell(xsf_file,gs,nz=nz,verbose=verbose)
  kpt = cell.get_abs_kpts([.0,.0,.0])  # gamma point calculation
  mf = scf.RHF(cell,exxdiv=exxdiv)     # exxdiv=None means no Ewald sum
  
  # tweak converger
  mf.conv_tol         = 1e-7
  mf.direct_scf_tol   = 1e-7
  mf.max_cycle        = 100
  mf.diis_start_cycle = 1

  # run or load
  mf.chkfile = chkfile_name # set restart file
  if os.path.isfile(chkfile_name):
    from pyscf import lib
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  return mf
# end def run_pyscf
