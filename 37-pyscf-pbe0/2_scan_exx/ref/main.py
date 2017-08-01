#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np
import pandas as pd
import h5py

def bfd_basis(nz=2):
  """ return a truncated BFD basis set for carbon
  nz = 2 double-zeta 
     = 3 triple-zeta etc. 
  """
  from pyscf.gto.basis import parse
  if nz == 2:
    mn_3d_start_idx=[12,46]
    basis={ # truncated BFD double-zeta basis set
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
  else:
    raise RuntimeError('nz must be one of 2 (double zeta) for now ')
  # end if
  return basis, mn_3d_start_idx
# end def bfd_basis

def build_cell(xsf_file,gs,nz=2,verbose=4):
  from pyscf.pbc import gto
  from nexus import Structure # use nexus.Structure to read xsf
  from pyscf_orbital_routines import atom_text # take ['C','C'],[[0,0,0],[0.5,0.5,0.5]] to make 'C 0 0 0\n C 0.5 0.5 0.5'

  # define Hamiltonian and atomic basis set
  pseudo = {'Mn':'bfd','O':'bfd'}
  basis,mn_3d_start_idx  = bfd_basis(nz=nz)

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
    unit    = 'anything_but_b_B_au_AU' # PySCF's default unit is angstrom
  )
  return cell,mn_3d_start_idx
# end def build_cell

def modify_1rdm(dm,mn_dorb_indices,ndorb=5):
  """ impose anti-ferromagnetic on UHF density matrix
  Inputs:
    dm: 3D numpy array of shape (nspin,nao,nao) - expect nspin=2 for up and down 
    mn_dorb_indices: a list of int, one for each Mn atom - index of the first d orbital
    ndorb: int, number of d orbitals - should be 5 right? 
  Output:
    new_dm: 3D numpy array, copied from dm, then modified """

  nspin,nao,nao1 = dm.shape
  assert nspin == 2
  new_dm = dm.copy()
  for ispin in [0,1]:

    # AFM order
    first_up = True
    if ispin == 0:
      first_up = False

    set_up = first_up
    for ibegin in mn_dorb_indices:
      # set all d electrons
      iend   = ibegin + ndorb
      if set_up:
        new_dm[ispin,ibegin:iend,ibegin:iend] = 1.0 * np.eye(ndorb)
      else:
        new_dm[ispin,ibegin:iend,ibegin:iend] = 0.0 * np.eye(ndorb)
      # end if
      set_up = not set_up
    # end for ibegin

  # end for ispin
  return new_dm
# end def modify_1rdm

def run_pyscf(xsf_file,gs,xc,nz=2,verbose=4,chkfile_name='bfd.h5',exxdiv=None):
  import os
  from pyscf.pbc import dft,scf

  cell,mn_3d_start_idx = build_cell(xsf_file,gs,nz=nz,verbose=verbose)
  kpt = cell.get_abs_kpts([.0,.0,.0])  # gamma point calculation
  mf = dft.UKS(cell)     # exxdiv=None means no Ewald sum, UKS: unexpected keyword argument 'exxdiv'
  
  # tweak converger
  mf.conv_tol         = 1e-7
  mf.direct_scf_tol   = 1e-7
  mf.max_cycle        = 100
  mf.diis_start_cycle = 1
  mf.xc = xc

  # run or load
  mf.chkfile = chkfile_name # set restart file
  if os.path.isfile(chkfile_name):
    from pyscf import lib
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  # end if
  return mf, mn_3d_start_idx
# end def run_pyscf

if __name__ == '__main__':
  import sys
  exx_fraction = float( sys.argv[1] )

  from pyscf_orbital_routines import uhf_multideterminant_spos, generate_pwscf_h5

  from datetime import datetime
  proj_id    = 'mno'
  verbosity  = 3
  grid_shape = [32]*3
  pseudos    = {'Mn':'Mn.BFD.xml','O':'O.BFD.xml'}
  ndet       = 1
  gvec_fname = 'gvectors.dat'
  eig_fname  = 'eigensystem.json'
  h5_fname   = 'pyscf2qmcpack.h5'
  vmc_inp    = 'vmc.xml'
  opt_inp    = 'opt.xml'
  dmc_inp    = 'dmc.xml'
  xsf_file   = '../../struct/mno.xsf'
  fockh5     = 'fock.h5'
  chkfile    = 'bfd.h5'
  nloop    = 5
  nwalker  = 4608 # 36*128
  xc = '%3.2f*HF+%3.2f*PBE,PBE' % (exx_fraction,1.-exx_fraction)
  skip_afm = True

  assert ndet == 1 # this is a single determinant run

  # run pyscf to obtain mean-field object
  print('running HF...')
  print(datetime.now())
  mf, mn_3d_start_idx = run_pyscf(xsf_file,grid_shape,xc,verbose=verbosity,chkfile_name=chkfile)
  print(datetime.now())
  print(' HF done')
  print(' SCF energy = %f'%mf.e_tot)

  # converge to AFM state
  if not skip_afm:
    dm = mf.make_rdm1()
    new_dm  = modify_1rdm(dm,mn_3d_start_idx)
    mf.max_cycle        = 200
    mf.kernel(new_dm)
    mp  = mf.mulliken_pop()
  # end if
  
  # save hcore and fock matrix for fcidump (next step phf)
  if not os.path.isfile(fockh5):
    sp.check_call(['cp',chkfile,fockh5])
    # evaluate hcore and fock matrix
    print('evaluating integrals...')
    print(datetime.now())
    hcore = mf.get_hcore()
    fock  = (hcore + mf.get_veff())
    print(datetime.now())
    print(' integrals done')
    with h5py.File(fockh5) as fh5:
      fh5['scf/hcore'] = hcore
      fh5['scf/fock']  = fock
    # end with
  # end if

  # save Kohn-Sham eigensystem to dataframe
  if not (os.path.isfile(gvec_fname) and os.path.isfile(eig_fname)):
    nup,ndn = mf.cell.nelec
    nspin,nao,nmo = mf.mo_coeff.shape # !!!! UHF, will also pass for kpoint, which is bad

    # initialize one determinant (Hatree-Fock)
    det_list = np.zeros([ndet,nmo,nmo,nspin])
    for ispin in range(nspin):
      det_list[0,:,:,ispin] = np.eye(nmo)
    # end for

    gvecs,eig_df = uhf_multideterminant_spos(det_list,{0:nup,1:ndn},
      mf.cell,mf.mo_coeff,mf.kpt)
    np.savetxt(gvec_fname,gvecs)
    eig_df.reset_index().to_json(eig_fname)
  # end if

  # write QMCPACK wavefunction file
  if not os.path.isfile(h5_fname):
    gvecs  = np.loadtxt(gvec_fname)
    eig_df = pd.read_json(eig_fname).set_index(['ikpt','ispin','istate']).sort_index()
    print('FFT orbitals ...')
    print(datetime.now())
    generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge={'Mn':10,'O':2},h5_fname=h5_fname)
    print(datetime.now())
    print(' FFT done')
  # end if

  # write QMCPACK input
  if True:#not (os.path.isfile(opt_inp) and os.path.isfile(dmc_inp)):
    cell    = mf.cell
    nup,ndn = cell.nelec
    fftgrid = 2*np.array(cell.gs)+1
    from input_xml import InputXml
    inp = InputXml()

    # UHF determinant - need to add Jastrow!
    wf_node = inp.uhf_slater(h5_fname,{'u':nup,'d':ndn}
      ,fftgrid=' '.join(fftgrid.astype(str)))

    # write vmc input
    vmc_node = inp.get_qmc_node(144)
    inp.write_qmcpack_input(vmc_inp,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_nodes=[vmc_node],proj_id=proj_id,pseudos=pseudos)

    # write opt input - need to add Jastrow!
    opt_node = inp.get_optimization_node(nloop)
    inp.write_qmcpack_input(opt_inp,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_nodes=[opt_node],proj_id=proj_id,pseudos=pseudos)

    # write dmc input - need to add Jastrow!
    dmc_nodes = inp.get_dmc_nodes(nwalker,nvmc_walkers=16,time_step_list=[0.008,0.004],correlation_time=0.5)
    inp.write_qmcpack_input(dmc_inp,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_nodes=dmc_nodes,proj_id=proj_id,pseudos=pseudos)
  # end if

# end __main__
