#!/usr/bin/env python
import os
import sys
sys.path.insert(0,'../2_eigensys')
from carbon import run_carbon
import numpy as np
from lxml import etree

if __name__ == '__main__':

  mf = run_carbon()
  nao = mf.mo_coeff.shape[0]

  ci_coeff = np.loadtxt('../1_ref/ci_coeff.dat').view(complex)
  detlist  = np.loadtxt('../1_ref/detlist.dat').view(complex)
  ndet = len(detlist)
  assert len(ci_coeff) == ndet
  assert detlist.shape[1] == nao*nao


  from input_xml import InputXml
  inp = InputXml()

  # build <simulationcell>
  sc_node = inp.simulationcell_from_cell(mf.cell)

  # build <particleset>
  import h5py
  fp = h5py.File('../2_eigensys/pwscf.pwscf.h5')
  pset_node = inp.particleset_from_hdf5(fp)

  # build <sposet>

  # build <determinantset>

  # build <hamiltonian>

  # assemble <qmcsystem>
  sys_node = etree.Element('qmcsystem')
  sys_children = [sc_node,pset_node]
  for child in sys_children:
    sys_node.append(child)
  # end for

  # write input
  from lxml import etree
  root = etree.Element('simulation')
  doc = etree.ElementTree(root)
  root.append(sys_node)
  doc.write('test.xml',pretty_print=True)
# end __main__
