#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np
import pandas as pd
import h5py

if __name__ == '__main__':

  import sys
  sys.path.insert(0,'scripts')
  from step1_run_pyscf import run_pyscf
  from pyscf_orbital_routines import multideterminant_orbitals, generate_pwscf_h5

  from datetime import datetime
  proj_id    = 'c2'
  verbosity  = 3
  grid_shape = [16]*3
  ndet       = 1
  gvec_fname = 'gvectors.dat'
  eig_fname  = 'eigensystem.json'
  h5_fname   = 'pyscf2qmcpack.h5'
  opt_inp    = 'opt.xml'
  dmc_inp    = 'dmc.xml'
  xsf_file   = '../../struct/c2.xsf'
  fockh5     = 'fock.h5'
  chkfile    = 'bfd.h5'

  assert ndet == 1 # this is a single determinant run

  # run pyscf to obtain mean-field object
  print('running HF...')
  print(datetime.now())
  mf = run_pyscf(xsf_file,grid_shape,verbose=verbosity,chkfile_name=chkfile)
  print(datetime.now())
  print(' HF done')

  ewald_yes = -10.2534689475175
  ewald_no  = -7.55768592061669
  ewald = ewald_yes-ewald_no
  print('ewald correction: %f' % ewald)

  checking = False # dump integrals in ascii file for old phfmol
  if checking:
    from pyscf import tools, ao2mo
    from functools import reduce

    mf.kernel()
    c = mf.mo_coeff
    nmo = c.shape[-1]
    h1e = reduce(np.dot, (c.T, mf.get_hcore(), c))
    eri = mf.with_df.ao2mo(c) # 2-electron integrals
    eri = ao2mo.restore('s8',eri,nmo)
    tools.fcidump.from_integrals('fcidump.dat',
      h1e, eri, nmo, mf.cell.nelectron,nuc=mf.cell.energy_nuc(), ms=0)
  # end if checking
  
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
    assert nup == ndn
    nfill   = nup # !!!! RHF
    nao,nmo = mf.mo_coeff.shape

    # initialize one determinant (Hatree-Fock)
    det_list = np.zeros([ndet,nmo,nmo])
    det_list[0,:,:] = np.eye(nmo)

    gvecs,eig_df = multideterminant_orbitals(det_list,nfill,mf.cell,mf.mo_coeff,mf.kpt)
    np.savetxt(gvec_fname,gvecs)
    eig_df.reset_index().to_json(eig_fname)
  # end if

  # write QMCPACK wavefunction file
  if not os.path.isfile(h5_fname):
    gvecs  = np.loadtxt(gvec_fname)
    eig_df = pd.read_json(eig_fname).set_index(['ikpt','ispin','istate']).sort_index()
    print('FFT orbitals ...')
    print(datetime.now())
    generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge={'C':2},h5_fname=h5_fname)
    print(datetime.now())
    print(' FFT done')
  # end if

  # write QMCPACK input
  if not (os.path.isfile(opt_inp) and os.path.isfile(dmc_inp)):
    cell    = mf.cell
    nup,ndn = cell.nelec
    assert nup == ndn # !!!! RHF
    fftgrid = 2*np.array(cell.gs)+1
    from input_xml import InputXml
    inp = InputXml()

    # RHF determinant - need to add Jastrow!
    wf_node = inp.rhf_slater(h5_fname,nup,fftgrid=' '.join(fftgrid.astype(str)))

    pseudos = {'C':'C.BFD.xml'}
    nloop = 5
    opt_node = inp.get_optimization_node(nloop)
    inp.write_qmcpack_input(opt_inp,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_nodes=[opt_node],proj_id=proj_id,pseudos=pseudos)

    # write dmc input - need to add Jastrow!
    nwalker  = 4608 # 36*128
    dmc_nodes = inp.get_dmc_nodes(nwalker,nvmc_walkers=16)
    inp.write_qmcpack_input(dmc_inp,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_nodes=dmc_nodes,proj_id=proj_id,pseudos=pseudos)
  # end if

# end __main__
