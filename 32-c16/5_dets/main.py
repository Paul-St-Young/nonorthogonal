#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np
import pandas as pd
import h5py

if __name__ == '__main__':

  import sys
  sys.path.insert(0,'../2_hf/scripts')
  from step1_run_pyscf import run_pyscf
  from pyscf_orbital_routines import multideterminant_orbitals, generate_pwscf_h5

  from datetime import datetime
  proj_id    = 'c2'
  verbosity  = 3
  grid_shape = [16]*3
  gvec_fname = 'gvectors.dat'
  eig_fname  = 'eigensystem.json'
  detlist_fname = 'det_list.dat'
  h5_fname   = 'pyscf2qmcpack.h5'
  vmc_inp    = 'vmc.xml'
  xsf_file   = 'c2-16.xsf'

  # run pyscf to obtain mean-field object
  print('running HF...')
  print(datetime.now())
  mf = run_pyscf(xsf_file,grid_shape,verbose=verbosity)
  print(datetime.now())
  print(' HF done')

  # read determinants
  det_list  = np.loadtxt(detlist_fname).view(complex)
  ndet,nmo2 = det_list.shape
  nup,ndn = mf.cell.nelec
  assert nup == ndn
  nfill   = nup # !!!! RHF
  nao,nmo = mf.mo_coeff.shape
  assert nmo2 == nmo*nmo

  # save Kohn-Sham eigensystem to dataframe
  if not (os.path.isfile(gvec_fname) and os.path.isfile(eig_fname)):

    gvecs,eig_df = multideterminant_orbitals(det_list.reshape(ndet,nmo,nmo),
      nfill,mf.cell,mf.mo_coeff,mf.kpt)
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

  def rhf_multidet_vmc(vmc_infile,ci_coeff):
    from input_xml import InputXml
    inp = InputXml()
    cell    = mf.cell
    nup,ndn = cell.nelec
    assert nup == ndn # !!!! RHF
    fftgrid = 2*np.array(cell.gs)+1
    wf_node = inp.rhf_slater(h5_fname,nup,fftgrid=' '.join(fftgrid.astype(str)))

    # !!!! hack to make 2 copies of <sposet>
    sb_node = wf_node.find('.//sposet_builder')
    ss_nodes = sb_node.findall('.//sposet')
    assert len(ss_nodes) == 1 # only 1 <sposet> for RHF
    ss_name = ss_nodes[0].get('name')
    from copy import deepcopy
    ss_node1 = deepcopy(ss_nodes[0])
    ss_node1.set('name',ss_name+'1')
    sb_node.append(ss_node1)
    ss_name1 = ss_node1.get('name')

    mdet_node = inp.multideterminant_from_ci(ci_coeff,nfill,nfill*ndet,real_coeff=True,spo_up=ss_name,spo_dn=ss_name1)

    # swap out <slaterdeterminant> for <multideterminant>
    ds_node = wf_node.find('.//determinantset')
    sd_node = ds_node.find('.//slaterdeterminant')
    ds_node.remove(sd_node)
    ds_node.append(mdet_node)

    # write vmc input - no Jastrow
    nwalker  = 144 # need a lot of walkers to fight variance from cusps
    vmc_node = inp.get_qmc_node(nwalker,checkpoint=-1)
    inp.write_qmcpack_input(vmc_infile,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_nodes=[vmc_node],proj_id=proj_id)

  # write QMCPACK inputs
  subdir = 'hii'
  if not os.path.isdir(subdir):
    sp.check_call(['mkdir',subdir])
  # end if
  for idet in range(ndet):
    detdir = os.path.join(subdir,'det%d'%idet)
    if not os.path.isdir(detdir):
      sp.check_call(['mkdir',detdir])
    # end if
    fake_ci_coeff = np.zeros(ndet)
    fake_ci_coeff[idet] = 1.0
    inp_loc = os.path.join(detdir,vmc_inp)
    rhf_multidet_vmc(inp_loc,fake_ci_coeff)
  # end for idet

# end __main__
