#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np
import pandas as pd
import h5py

def build_cell(gs,basis,pseudo={'Mn':'bfd','O':'bfd'},verbose=4):
  from pyscf.pbc import gto
  from pyscf_orbital_routines import atom_text # take ['C','C'],[[0,0,0],[0.5,0.5,0.5]] to make 'C 0 0 0\n C 0.5 0.5 0.5'

  # !!!! hard code MnO 4-atom unit cell
  axes = np.array(
  [[ 4.42764957,  2.21382478,  2.21382478],
   [ 2.21382478,  4.42764957,  2.21382478],
   [ 2.21382478,  2.21382478,  4.42764957]]
  )
  pos  = np.array(
  [[ 0.        ,  0.        ,  0.        ],
   [ 4.42764957,  4.42764957,  4.42764957],
   [ 2.21382478,  2.21382478,  2.21382478],
   [ 6.64147435,  6.64147435,  6.64148026]]
  )
  elem = ['Mn','Mn','O','O']

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
  return cell
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

def run_pyscf(gs,basis,verbose=4,chkfile_name='bfd.h5',exxdiv=None):
  import os
  from pyscf.pbc import dft,scf

  cell = build_cell(gs,basis,verbose=verbose)
  kpt  = cell.get_abs_kpts([.0,.0,.0])  # gamma point calculation
  mf   = scf.UHF(cell,exxdiv=exxdiv)    # exxdiv=None means no Ewald sum, UKS: unexpected keyword argument 'exxdiv'
  
  # tweak converger
  mf.conv_tol         = 1e-7
  mf.direct_scf_tol   = 1e-7
  mf.max_cycle        = 50 # probably not enough to converge, intended to be reran using a different guess (e.g. AFM DM)
  mf.diis_start_cycle = 1
  mf.max_memory       = 128000 # 128 GB

  # run or load
  mf.chkfile = chkfile_name # set restart file
  if os.path.isfile(chkfile_name):
    from pyscf import lib
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  return mf
# end def run_pyscf

if __name__ == '__main__':

  from pyscf_orbital_routines import uhf_multideterminant_spos, generate_pwscf_h5
  from datetime import datetime

  ndet       = 1 # single-deterimant calculation

  # PySCF inputs
  nz         = 2 # double-zeta BFD basis
  mygs       = 40
  gs         = [mygs]*3
  chkfile    = 'gs%s.h5' % mygs
  verbosity  = 4

  # QMCPACK converter inputes
  proj_id    = 'mno'
  pseudos    = {'Mn':'Mn.BFD.xml','O':'O.BFD.xml'}
  fockh5     = 'fock.h5' # hdf5 to store integrals (scf/hcore & scf.fock)
  gvec_fname = 'gvectors.dat'
  eig_fname  = 'eigensystem.json'
  h5_fname   = 'pyscf2qmcpack.h5'
  vmc_inp    = 'vmc.xml'
  opt_inp    = 'opt.xml'
  dmc_inp    = 'dmc.xml'
  nloop    = 5
  nwalker  = 4608 # 36*128

  assert ndet == 1 # this is a single determinant run

  # get PySCF basis
  import sys
  sys.path.insert(0,'../../../utils')
  from bfd_basis import mno_basis
  basis, mn_3d_idx = mno_basis(nz=nz)
  # locations of the 3d orbitals are needed to set AFM density matrix

  # run pyscf to obtain mean-field object (probably not yet converged)
  print('  running HF...')
  print(datetime.now())
  mf = run_pyscf(gs,basis,verbose=verbosity,chkfile_name=chkfile)
  print(datetime.now())
  print(' HF done')
  print(' SCF energy = %f'%mf.e_tot)

  print('restarting with modified (AFM) density matrix')
  dm = mf.make_rdm1()
  new_dm  = modify_1rdm(dm,mn_3d_idx)
  mf.max_cycle = 200
  mf.kernel(new_dm)
  mp  = mf.mulliken_pop()

  assert mf.converged, " SCF cycled is not converged .... "
  
  # save hcore and fock matrix for fcidump (next step phf)
  if not os.path.isfile(fockh5):
    sp.check_call(['cp',chkfile,fockh5])
    # evaluate hcore and fock matrix
    print('  evaluating integrals...')
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
    print('  FFT orbitals ...')
    print(datetime.now())
    generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge={'Mn':10,'O':2},h5_fname=h5_fname)
    print(datetime.now())
    print(' FFT done')
  # end if

  # write QMCPACK input
  if not (os.path.isfile(opt_inp) and os.path.isfile(dmc_inp)):
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
