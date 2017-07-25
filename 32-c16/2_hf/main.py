#!/usr/bin/env python
import os
import numpy as np
import pandas as pd

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
  xsf_file   = 'c2-16.xsf'

  assert ndet == 1 # this is a single determinant run

  # run pyscf to obtain mean-field object
  print('running HF...')
  print(datetime.now())
  mf = run_pyscf(xsf_file,grid_shape,verbose=verbosity)
  print(datetime.now())
  print(' HF done')

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
