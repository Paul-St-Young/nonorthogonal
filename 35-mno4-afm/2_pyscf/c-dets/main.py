#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np
import pandas as pd
import h5py

if __name__ == '__main__':

  from pyscf.pbc import scf
  from datetime import datetime
  from pyscf_orbital_routines import uhf_multideterminant_spos, generate_pwscf_h5

  proj_id    = 'mno'
  pseudos    = {'Mn':'Mn.BFD.xml','O':'O.BFD.xml'}
  pseudized_charges = {'Mn':10,'O':2} # number of electrons removed by the pseudopotential
  gvec_fname = 'gvectors.dat'
  eig_fname  = 'eigensystem.json'
  detlist_fname = '../b-phf/det_list.dat'
  vmc_inp    = 'vmc.xml'
  chkfile    = '../a-hf/bfd.h5'
  LINDEP_CUTOFF = 1e-9
  fftgrid = np.array([100,100,100])

  h5_loc  = 'pyscf2qmcpack.h5'

  # construct mf from chkfile
  cell,scf_rec = scf.chkfile.load_scf(chkfile)
  mf = scf.RHF(cell)
  mf.__dict__.update(scf_rec)

  # read determinants
  det_list  = np.loadtxt(detlist_fname).view(complex)
  ndet,nmo2 = det_list.shape
  nup,ndn = mf.cell.nelec
  assert nup == ndn

  mo_coeff = mf.mo_coeff
  nspin,nao,nmo = mo_coeff.shape

  useX = True
  if useX: # use X matrix
    from my_qmctools.new_integrals_from_chkfile import getOrthoAORotation
    X,nmo_per_kpt = getOrthoAORotation(cell,mf.kpt,LINDEP_CUTOFF)

    for ispin in range(nspin):
      mo_coeff[ispin] = X.copy()
    # end for ispin
    print('useful X matrix columns = %d' % nmo_per_kpt)
  # end if

  assert nmo2 == nmo*nmo*nspin, '%d entries in determinant, %d orbitals %d spin -> expect %d' % (nmo2,nmo,nspin,nmo*nmo*nspin)

  # save Kohn-Sham eigensystem to dataframe
  if not (os.path.isfile(gvec_fname) and os.path.isfile(eig_fname)):
    gvecs,eig_df = uhf_multideterminant_spos(det_list.reshape(ndet,nmo,nmo,nspin),
      {0:nup,1:ndn},mf.cell,mo_coeff,mf.kpt)
    np.savetxt(gvec_fname,gvecs)
    eig_df.reset_index().to_json(eig_fname)
  # end if

  # write QMCPACK wavefunction file
  if not os.path.isfile(h5_loc):
    gvecs  = np.loadtxt(gvec_fname)
    eig_df = pd.read_json(eig_fname).set_index(['ikpt','ispin','istate']).sort_index()
    print('FFT orbitals ...')
    print(datetime.now())
    generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge=pseudized_charges,h5_fname=h5_loc)
    print(datetime.now())
    print(' FFT done')
  # end if

  # UHF multi-determinant - need to add Jastrow!
  from input_xml import InputXml
  inp = InputXml()
  nwalker  = 144 # need a lot of walkers to fight variance from cusps
  vmc_node = inp.get_qmc_node(nwalker,checkpoint=-1)

  #   write QMCPACK inputs
  subdir = 'hii'
  if not os.path.isdir(subdir):
    sp.check_call(['mkdir',subdir])
  # end if
  for idet in range(ndet):
    detdir = os.path.join(subdir,'det%d'%idet)
    if not os.path.isdir(detdir):
      sp.check_call(['mkdir',detdir])
    # end if
    fake_ci_coeff = np.zeros(idet+1)
    fake_ci_coeff[idet] = 1.0
    inp_loc  = os.path.join(detdir,vmc_inp)
    h5_href  = os.path.relpath(h5_loc, os.path.dirname(inp_loc) )
    wf_node  = inp.uhf_multidet_qmc(fake_ci_coeff,nup,ndn,fftgrid,h5_href=h5_href,real_coeff=True)
    inp.write_qmcpack_input(inp_loc,mf.cell,h5_loc,nup,ndn,wf_node=wf_node,pseudos=pseudos,qmc_nodes=[vmc_node],proj_id=proj_id)
  # end for idet

# end __main__
