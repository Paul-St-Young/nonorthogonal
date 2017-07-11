import numpy as np
import h5py
from pyscf.pbc import gto, scf

if __name__ == '__main__':
  
  alat = 2.34301125

  axes_text = '''
   {alat:12.10f}  {half:12.10f}  {half:12.10f}
   {half:12.10f}  {alat:12.10f}  {half:12.10f}
   {half:12.10f}  {half:12.10f}  {alat:12.10f}
  '''.format(alat=alat,half=alat/2.)
  atom_text = '''
   Mn      0.00000      0.00000      0.00000
   Mn      {alat:12.10f} {alat:12.10f} {alat:12.10f}
   O       {half:12.10f} {half:12.10f} {half:12.10f}
   O       {half3:12.10f} {half3:12.10f} {half3:10.6}
  '''.format(alat=alat,half=alat/2.,half3=3.*alat/2.)

  cell = gto.M(
    verbose = 4,
    gs      = [16]*3,
    a       = axes_text,
    atom    = atom_text,
    pseudo  = {'Mn':'bfd','O':'bfd'},
  )

  from nexus import Structure
  assert cell.unit == 'angstrom'
  struct = Structure(axes = cell.lattice_vectors(),elem=['Mn','Mn','O','O'],pos=cell.atom_coords(),units='A')
  struct.write('mno.xsf')

# end __main__
