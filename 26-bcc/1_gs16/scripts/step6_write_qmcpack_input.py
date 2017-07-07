def step6_write_qmcpack_input(inp_name,cell,wf_h5_fname,nup,ndn,proj_id='c2'):
  from input_xml import InputXml
  inp = InputXml()

  # build <project>
  from lxml import etree
  proj_node = etree.Element('project',{'id':'c2','series':'0'})

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
    sys_node.append(child)
  # end for

  # write <qmc> block for a quick VMC
  nblock = 400
  nstep  = 10
  time_step = 2.0
  nwalker   = 16
  blocks_node = etree.Element('parameter',{'name':'blocks'})
  blocks_node.text = str(nblock)
  steps_node = etree.Element('parameter',{'name':'steps'})
  steps_node.text = str(nstep)
  ts_node = etree.Element('parameter',{'name':'timestep'})
  ts_node.text = str(time_step)
  walker_node = etree.Element('parameter',{'name':'walkers'})
  walker_node.text = str(nwalker)
  vmc_children = [blocks_node,steps_node,ts_node,walker_node]
  vmc_node = etree.Element('qmc',{'method':'vmc','move':'pbyp'})
  for child in vmc_children:
    vmc_node.append(child)
  # end for

  # write input
  from lxml import etree
  root = etree.Element('simulation')
  doc = etree.ElementTree(root)
  root.append(sys_node)
  root.append(vmc_node)
  doc.write(inp_name,pretty_print=True)
# end def step6_write_qmcpack_input
