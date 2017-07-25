#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import h5py

def write_qmcpack_input(inp_name,cell,wf_h5_fname,nup,ndn,wf_node=None,qmc_nodes=[],proj_id='qmc'):
  from input_xml import InputXml
  inp = InputXml()

  # build <project>
  from lxml import etree
  proj_node = etree.Element('project',{'id':proj_id,'series':'0'})

  # build <simulationcell>
  sc_node   = inp.simulationcell_from_cell(cell)

  # build <particleset>
  elec_pset_node= inp.ud_electrons(nup,ndn)
  import h5py
  fp = h5py.File(wf_h5_fname)
  ion_pset_node = inp.particleset_from_hdf5(fp)

  # build <wavefunction>
  #  in another file

  # build <hamiltonian>
  ii_node = etree.Element('constant',{'type':'coulomb','name':'IonIon'
    ,'source':ion_pset_node.get('name'),'target':ion_pset_node.get('name')})
  ee_node = etree.Element('pairpot',{'type':'coulomb','name':'ElecElec'
    ,'source':elec_pset_node.get('name'),'target':elec_pset_node.get('name')})

  # !!!! hard-code electron-ion pseudized interaction
  ei_node = etree.Element('pairpot',{'type':'pseudo','name':'PseudoPot'
    ,'source':ion_pset_node.get('name'),'target':elec_pset_node.get('name')
    ,'wavefunction':'psi0','format':'xml'})
  pseudo_node = etree.Element('pseudo',{'elementType':'C','href':'C.BFD.xml'})
  ei_node.append(pseudo_node)

  ham_children = [ii_node,ee_node,ei_node]
  ham_node = etree.Element('hamiltonian',{'name':'h0','type':'generic','target':elec_pset_node.get('name')})
  for child in ham_children:
    ham_node.append(child)
  # end for

  # assemble <qmcsystem>
  sys_node = etree.Element('qmcsystem')
  sys_children = [proj_node,sc_node,elec_pset_node,ion_pset_node,ham_node]

  for child in sys_children:
    # if give, insert <wavefunction> before <hamiltonian> 
    if (child.tag == 'hamiltonian') and (wf_node is not None):
      sys_node.append(wf_node)
    # end if
    sys_node.append(child)
  # end for

  # write input
  from lxml import etree
  root = etree.Element('simulation')
  doc = etree.ElementTree(root)
  root.append(sys_node)

  # take <qmc> block from else where
  if len(qmc_nodes) > 0:
    for qmc_node in qmc_nodes:
      root.append(qmc_node)
    # end for
  # end if

  doc.write(inp_name,pretty_print=True)
# end def write_qmcpack_input

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
    wf_node = inp.rhf_slater(h5_fname,nup,fftgrid=' '.join(fftgrid.astype(str)))

    # write optimization input - need to add Jastrow!
    nloop = 5
    opt_node = inp.get_optimization_node(nloop)
    write_qmcpack_input(opt_inp,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_nodes=[opt_node],proj_id=proj_id)

    # write dmc input - need to add Jastrow!
    nwalker  = 4608 # 36*128
    dmc_nodes = inp.get_dmc_nodes(nwalker,nvmc_walkers=16)
    write_qmcpack_input(dmc_inp,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_nodes=dmc_nodes,proj_id=proj_id)
  # end if

# end __main__
