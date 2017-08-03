#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np
from lxml import etree

def main():
  inp_loc = './scan/ldau/lda-0-ecut320/u0.50/dmc/lda-0-ecut320-dmc.in.xml'
  cont_iseries = 2
  is_text = str(cont_iseries).zfill(3)
  cont_postfix = '.s%s.config.h5' % is_text  # walkers
  rand_postfix = '.s%s.random.xml' % is_text # random number generator state
  qmc_postfix  = '.s%s.qmc.xml' % is_text    # branch engine state
  sub_postfix = '.sbatch.in'
  ntimes = 16 # multiply the number of blocks
  pseudos = ['Mn.BFD.xml','O.BFD.xml']

  from input_xml import InputXml
  inp = InputXml()
  inp.read(inp_loc)

  # get project id
  proj_node = inp.find('.//project')
  prefix = proj_node.get('id')
  proj_node.set('series','1') # skipping VMC

  # add estimators
  ham = inp.find('.//hamiltonian')
  gr_num_bin = 50
  grid_shape = np.array([20,20,20])
  for est in [inp.sk(),inp.gr(gr_num_bin),inp.spindensity(grid_shape),inp.sksp(),inp.localmoment()]:
    ham.append( est )
  # end for

  # read previous walkers
  restart_node = etree.Element('mcwalkerset',{'fileroot':prefix+'.s'+is_text,'node':'-1','nprocs':'576','version':'3 1','collected':'yes'})
  first_qmc_node = inp.find('.//qmc')
  sim = first_qmc_node.getparent()
  idx = sim.index(first_qmc_node)
  sim.insert(idx,restart_node)

  # edit qmc section
  for qmc_node in inp.findall('.//qmc'):
    method = qmc_node.get('method')
    if method == 'vmc':
      sim.remove(qmc_node)
    else:
      nblock_node = qmc_node.find('.//parameter[@name="blocks"]')
      nblock = int( nblock_node.text )
      nblock_node.text = str(nblock*ntimes)
    # end if
  # end for

  # find wavefunction file
  inp_dir= os.path.dirname(inp_loc)
  inp2h5 = inp.find('.//sposet_builder').get('href')
  h5_loc = os.path.abspath( os.path.join(inp_dir,inp2h5) )

  # make new run folder
  new_dir = inp_dir.replace('./scan/','./cont_scan/')
  sp.check_call(['mkdir','-p',new_dir])
  h5_rel = os.path.relpath(h5_loc,new_dir)
  inp.find('.//sposet_builder').set('href',h5_rel) # edit h5 reference

  # send restart file
  cont_file = os.path.join(inp_dir,prefix+cont_postfix)
  sp.check_call(['cp',cont_file,new_dir])
  # send random number file
  rand_file = os.path.join(inp_dir,prefix+rand_postfix)
  sp.check_call(['cp',rand_file,new_dir])
  # send qmc file
  qmc_file = os.path.join(inp_dir,prefix+qmc_postfix)
  sp.check_call(['cp',qmc_file,new_dir])

  # send submission script
  sub_file = os.path.join(inp_dir,prefix+sub_postfix)
  new_sub  = os.path.join(new_dir,os.path.basename(sub_file))
  sp.check_call(['cp',sub_file,new_sub])
  # edit submission script
  with open(new_sub,'r') as f:
    sub_text = f.read()
  # end if
  sub_text.replace('-t 01:00:00','-t 04:00:00')
  with open(new_sub,'w') as f:
    f.write(sub_text)
  # end if

  # send pseudopotentials
  for psp in pseudos:
    sp.check_call(['cp',os.path.join(inp_dir,psp),os.path.join(new_dir,psp)])
  # end for

  # send new input
  new_infile = os.path.join( new_dir, os.path.basename(inp_loc) )
  inp.write(new_infile)

# end def main

if __name__ == '__main__':
  main()
