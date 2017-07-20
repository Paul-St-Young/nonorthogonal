#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import h5py

def modify_1rdm(dm,mn_dorb_indices,ndorb=5):
  """ impose anti-ferromagnetic order on UHF density matrix
  Inputs:
    dm: 3D numpy array of shape (nspin,nao,nao) - expect nspin=2 for up and down 
    mn_dorb_indices: a list of int, one for each Mn atom - index of the first d orbital
    ndorb: int, number of d orbitals - should be 5 for 3d right? 
  Output:
    new_dm: 3D numpy array, copied from dm, then modified """

  nspin,nao,nao1 = dm.shape
  assert nspin == 2
  new_dm = dm.copy()
  for ispin in [0,1]:

    # up and down spins cannot occupy the same orbital
    first_up = True
    if ispin != 0:
      first_up = False

    # impose AFM order
    occupy = first_up
    for ibegin in mn_dorb_indices: # loop through each Mn atom
      # set all 3d electrons
      iend   = ibegin + ndorb
      if occupy: # occupy orbitals
        new_dm[ispin,ibegin:iend,ibegin:iend] = 1.0 * np.eye(ndorb) 
      else: # vacate orbitals
        new_dm[ispin,ibegin:iend,ibegin:iend] = 0.0 * np.eye(ndorb) 
      # end if
      occupy = not occupy # AFM order i.e. neighboring atoms anti-align
    # end for ibegin

  # end for ispin
  return new_dm
# end def modify_1rdm

def write_uhf_multideterminant_spo(det_list,ndet,nfill,mf):
  import pandas as pd
  from pyscf_orbital_routines import multideterminant_orbitals

  # ndet0: number of available determinants
  ndet0,nspin,nmo,nmo1 = det_list.shape

  # split into up and down determinants to reuse RHF routines
  det_list_up = np.zeros([ndet,nmo,nmo],dtype=complex)
  det_list_dn = np.zeros([ndet,nmo,nmo],dtype=complex)
  for ispin in range(nspin): # up and down
    for idet in range(ndet):
      if ispin == 0:
        det_list_up[idet,:,:] = det_list[idet,ispin,:,:]
      elif ispin == 1:
        det_list_dn[idet,:,:] = det_list[idet,ispin,:,:]
      else:
        raise RuntimeError('cannot handle ispin=%d'%ispin)
      # end if
    # end for idet
  # end for ispin

  gvecs_up, eig_df_up = multideterminant_orbitals(det_list_up,nfill,mf.cell,mf.mo_coeff[0],mf.kpt,ikpt=0,ispin=0)
  gvecs_dn, eig_df_dn = multideterminant_orbitals(det_list_dn,nfill,mf.cell,mf.mo_coeff[1],mf.kpt,ikpt=0,ispin=1)
  assert np.allclose(gvecs_up,gvecs_dn)
  eig_df = pd.concat([eig_df_up,eig_df_dn])
  return gvecs_up,eig_df
# end def write_uhf_multideterminant_spo

def write_qmcpack_input(inp_name,cell,wf_h5_fname,nup,ndn,wf_node=None,qmc_node=None,proj_id='mno'):
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
  pseudo_node = etree.Element('pseudo',{'elementType':'Mn','href':'Mn.BFD.xml'})
  ei_node.append(pseudo_node)
  pseudo_node1 = etree.Element('pseudo',{'elementType':'O','href':'O.BFD.xml'})
  ei_node.append(pseudo_node1)
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

  # take <qmc> block from else where

  # write input
  from lxml import etree
  root = etree.Element('simulation')
  doc = etree.ElementTree(root)
  root.append(sys_node)
  if qmc_node is not None:
    root.append(qmc_node)
  doc.write(inp_name,pretty_print=True)
# end def write_qmcpack_input

if __name__ == '__main__':

  import sys
  sys.path.insert(0,'scripts')
  from step1_run_pyscf import run_pyscf

  from datetime import datetime
  skip_afm   = False
  verbosity  = 4
  grid_shape = [25]*3
  ndet       = 1
  gvec_fname = 'gvectors.dat'
  eig_fname  = 'eigensystem.json'
  h5_fname   = 'pyscf2qmcpack.h5'
  xml_fname  = 'qmc.xml'
  assert ndet == 1 # this is a single determinant run

  # run pyscf to obtain mean-field object
  print('running HF...')
  print(datetime.now())
  mf = run_pyscf(grid_shape,verbose=verbosity)
  print(datetime.now())
  print(' HF done')
  if not skip_afm:

    # get current 1RDM
    dm = mf.make_rdm1()
    np.savetxt('1rdm_up.dat',dm[0])
    np.savetxt('1rdm_dn.dat',dm[1])

    # enforce anti-ferromagnetic order
    new_dm = modify_1rdm(dm,[12,46])
    np.savetxt('new_1rdm_up.dat',new_dm[0])
    np.savetxt('new_1rdm_dn.dat',new_dm[1])

    # rerun scf
    print('rerunning HF...')
    print(datetime.now())
    mf.kernel(new_dm)
    print(datetime.now())
    print(' HF done')

    ## check state!
    print("Mulliken population:")
    mp  = mf.mulliken_pop()
    #labels = mf.cell.ao_labels(fmt=None)
    #with open('basis.txt','w') as fout:
    #  fout.write( ['%s  %s\n' % (label[2],label[3]) for label in labels] )
    ## end with
  # end if skip_afm

  # save Kohn-Sham eigensystem to dataframe
  if not (os.path.isfile(gvec_fname) and os.path.isfile(eig_fname)):
    nspin,nao,nmo = mf.mo_coeff.shape
    nfill = mf.cell.tot_electrons()

    # initialize one determinant (Hatree-Fock)
    det_list = np.zeros([ndet,nspin,nmo,nmo])
    for ispin in range(nspin):
      det_list[0,ispin,:,:] = np.eye(nmo)
    # end for

    gvecs,eig_df = write_uhf_multideterminant_spo(det_list,ndet,nfill,mf)
    np.savetxt(gvec_fname,gvecs)
    eig_df.reset_index().to_json(eig_fname)
  # end if

  # write QMCPACK wavefunction file
  if not os.path.isfile(h5_fname):
    gvecs  = np.loadtxt(gvec_fname)
    eig_df = pd.read_json(eig_fname).set_index(['ikpt','ispin','istate']).sort_index()
    from pyscf_orbital_routines import generate_pwscf_h5
    print('FFT orbitals ...')
    print(datetime.now())
    generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge={'Mn':10,'O':2},h5_fname=h5_fname)
    print(datetime.now())
    print(' FFT done')
  # end if

  # write QMCPACK input
  if not os.path.isfile(xml_fname):
    # mf is not really needed, can get cell using build_cell
    #from step1_run_pyscf import build_cell
    #cell    = build_cell(grid_shape,verbose=verbosity)
    cell    = mf.cell
    nup,ndn = cell.nelec
    fftgrid = 2*np.array(cell.gs)+1
    from input_xml import InputXml
    inp = InputXml()
    wf_node = inp.uhf_slater(h5_fname,{'u':nup,'d':ndn},fftgrid=' '.join(fftgrid.astype(str)))
    nloop = 5
    opt_node = inp.get_optimization_node(nloop)
    write_qmcpack_input(xml_fname,cell,h5_fname,nup,ndn,wf_node=wf_node,qmc_node=opt_node)
  # end if

# end __main__
