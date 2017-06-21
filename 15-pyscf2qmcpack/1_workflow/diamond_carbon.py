#!/usr/bin/env python
import numpy as np
import sys
sys.path.insert(0,'scripts')

def build_carbon_cell(verbose=4):
  from pyscf.pbc import gto
  import pyscf.gto.basis as bas
  basis={ # truncated BFD double-zeta basis set from Lucas K. Wagner
  'C':bas.parse('''
  C s
    0.205100 0.397529 
    0.409924 0.380369 
    0.819297 0.180113 
    1.637494 -0.033512 
    3.272791 -0.121499 
    6.541187 0.015176 
    13.073594 -0.000705 
   C p
    0.234064 0.302667 
    0.468003 0.289868 
    0.935757 0.210979 
    1.871016 0.112024 
    3.741035 0.054425 
    7.480076 0.021931 
   C d
    0.329486 1.000000 
   C f
    0.773485 1.000000 
   C s
  0.2 1.0
   C s
  0.6 1.0
   C p
  0.2 1.0
   C p
  0.6 1.0
  '''),
  }
  pseudo = {'C':'bfd'}
  alat0  = 3.6 # angstrom
  gs     = [4,4,4]

  # convert elem,pos to text representation
  def atom_text(elem,pos):
    assert len(elem) == len(pos)
    lines = []
    for iatom in range(len(elem)):
        mypos = pos[iatom]
        line = '%5s  %10.6f  %10.6f  %10.6f' % (elem[iatom],mypos[0],mypos[1],mypos[2])
        lines.append(line)
    atext = ';\n'.join(lines)
    return atext
  # end def
  axes  = (np.ones((3,3))-np.eye(3))*alat0/2.0
  elem  = ['C','C']
  pos   = np.array([[0,0,0],[0.25,0.25,0.25]])*alat0
  atoms = atom_text(elem,pos)
  cell = gto.M(verbose=verbose,a=axes,
    gs=gs,atom=atoms,basis=basis,pseudo=pseudo)
  return cell
# end def 

def run_carbon(verbose=4,chkfile_name='vdz.h5'):
  import os
  from pyscf.pbc.scf import RHF
  cell = build_carbon_cell(verbose=verbose)
  mf = RHF(cell)
  mf.max_cycle=50
  mf.chkfile=chkfile_name
  mf.conv_tol=1e-07
  mf.diis_start_cycle=1

  # run or load
  if os.path.isfile(chkfile_name):
    from pyscf import lib
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  return mf
# end def run_carbon

if __name__ == '__main__':

  # run pyscf and extract Kohn-Sham eigensystem
  # ================================================
  mean_field_object = run_carbon(verbose=3)

  from pyscf_orbital_routines import save_eigensystem
  gvecs, eig_df = save_eigensystem(mean_field_object,save=False)

  # generate wave function file
  # ================================================
  h5_fname = 'pyscf2pwscf.h5'
  from pyscf_orbital_routines import generate_pwscf_h5
  generate_pwscf_h5(mean_field_object.cell,gvecs,eig_df
    ,pseudized_charge={'C':2},h5_fname=h5_fname)

  # generate QMCPACK input file
  # ================================================
  xml_fname = 'test.xml'
  import h5py
  from lxml import etree
  from input_xml import InputXml
  h5_handle = h5py.File(h5_fname,'r')
  inp = InputXml()
  # build <simulationcell>
  sc_node = inp.simulationcell_from_cell(mean_field_object.cell)
  # build <particleset>
  pset_node = inp.particleset_from_hdf5(h5_handle)
  # assemble <qmcsystem>
  sys_node = etree.Element('qmcsystem')
  sys_children = [sc_node,pset_node]
  for child in sys_children:
    sys_node.append(child)
  # end for

  # write input
  root = etree.Element('simulation')
  doc = etree.ElementTree(root)
  root.append(sys_node)
  doc.write(xml_fname,pretty_print=True)
# end __main__
