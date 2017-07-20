import numpy as np

def build_cell(gs,verbose=4):
  from pyscf.pbc import gto
  from pyscf.gto.basis import parse

  alat = 4.427649569275 # angstrom
  pseudo = {'Mn':'bfd','O':'bfd'}
  basis={
  'Mn':parse('''
  Mn s
    23.6450680 -0.014659 
    13.4620290 0.223661 
    8.21262900 -0.564535 
    1.85994300 0.554600 
    0.88975400 0.534536 
    0.42256000 0.142766 
   Mn p
    17.8033170 0.003171 
    8.83499500 -0.080060 
    4.16134000 0.185735 
    2.12919600 0.375367 
    1.06371300 0.388624 
    0.51014100 0.192450 
    0.21130200 0.026146 
   Mn d
    9.25373600 0.076285 
    4.36301600 0.247722 
    1.83482100 0.330050 
    0.75611800 0.334518 
    0.29609200 0.263545 
   Mn f
    1.0998003483 1.00000 
   Mn s
  0.2 1.0
   Mn s
  0.6000000000000001 1.0
   Mn p
  0.2 1.0
   Mn p
  0.6000000000000001 1.0
   Mn d
  0.2 1.0
   Mn d
  0.6000000000000001 1.0
  '''),
  'O':parse('''
  O s
    0.268022 0.304848 
    0.573098 0.453752 
    1.225429 0.295926 
    2.620277 0.019567 
    5.602818 -0.128627 
    11.980245 0.012024 
    25.616801 0.000407 
    54.775216 -0.000076 
   O p
    0.333673 0.255999 
    0.666627 0.281879 
    1.331816 0.242835 
    2.660761 0.161134 
    5.315785 0.082308 
    10.620108 0.039899 
    21.217318 0.004679 
   O d
    0.669340 1.000000 
   O f
    1.423104 1.000000 
   O s
  0.2 1.0
   O s
  0.6000000000000001 1.0
   O p
  0.2 1.0
   O p
  0.6000000000000001 1.0
  '''),
  }

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
    verbose = verbose,
    gs      = gs,
    a       = axes_text,
    atom    = atom_text,
    basis   = basis,
    pseudo  = pseudo,
    unit    = 'angstrom'
  )
  return cell
# end def build_cell

def run_pyscf(gs,verbose=4,chkfile_name='bfd.h5'):
  import os
  from pyscf.pbc import scf

  cell = build_cell(gs,verbose)
  kpt = cell.get_abs_kpts([.0,.0,.0])  # gamma point calculation
  mf = scf.UHF(cell,exxdiv=None)       # exxdiv=None means no Ewald sum
  
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
