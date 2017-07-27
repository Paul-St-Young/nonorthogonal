#!/usr/bin/env python
# launch qmc calculations using template input files written by main.py
import os
import subprocess as sp
from lxml import etree
from input_xml import InputXml

def submission_template():
  ftext = '''#!/bin/bash
#SBATCH --job-name {name:s}
#SBATCH --tasks-per-node={npp:d}
#SBATCH --cpus-per-task=1
#SBATCH -N {nodes:d}
#SBATCH -p {queue:s}
#SBATCH -t {minutes:d}

BIN={bin_loc:s}

export OMP_NUM_THREADS={threads:d}
date
srun -n {nproc:d} $BIN {inp:s} > out 2> err
date '''
  return ftext
# end def submission_template

def run_inp_in_mydir(inp,mydir,sub_sbatch,pseudo_dir,pseudo_files,submit=True,name='qmc'):
  # make dir
  if not os.path.isdir(mydir):
    sp.check_call(['mkdir',mydir])
  # end if

  # send input
  sp.check_call(['cp',inp,mydir])
  inp_loc = os.path.join(mydir,inp)
  doc = etree.parse(inp_loc)
  root = doc.getroot()
  h5_node = root.find('.//sposet_builder')
  wf_h5_fname = h5_node.get('href')
  h5_node.set('href',os.path.join('../',wf_h5_fname))
  doc.write(inp_loc)

  # write and send submission script
  sub_text = submission_template()
  nodes = 1
  npp   = 36
  nproc = nodes*npp
  sub_content = sub_text.format(
    name = name,
    nodes= nodes,
    npp  = npp,
    nproc= nproc,
    queue= 'pbatch',
    minutes = 240,
    bin_loc = '/g/g91/yang41/soft/quartz_qmcpack/build/bin/qmcpack',
    threads = 1,
    inp     = inp
  )
  sub_loc = os.path.join(mydir,sub_sbatch)
  with open(sub_loc,'w') as f:
    f.write(sub_content)
  # end with

  # send pseudopotentials
  for psp in pseudo_files:
    psp_src = os.path.join(pseudo_dir,psp)
    cmd = ' '.join(['cp',psp_src,mydir])
    #sp.check_call(cmd)
    os.system(cmd)
  # end for

  out = ''
  err = ''
  if submit:
    cmd  = 'cd {subdir:s};sbatch {sub_file:s}'.format(subdir=mydir,sub_file=sub_sbatch)
    proc = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
    out,err = proc.communicate()
  # end if submit

  return out,err
# end run_inp_in_mydir

if __name__ == '__main__':

  run_opt = True
  run_dmc = False
  opt_wf_fname = 'jopt/c2.s004.opt.xml'

  pseudo_dir   = '~/pseudo'
  pseudo_files = ['C.BFD.xml']
  sub_sbatch   = 'qmc.sbatch' # name of submision file, written from scratch

  opt_dir = 'jopt'
  opt_inp = 'opt.xml'
  jas_opt_inp = 'jopt.xml'

  dmc_dir = 'dmc'
  dmc_inp = 'dmc.xml'
  jas_dmc_inp = 'jdmc.xml'

  out_dump = 'loc_launch.txt'
  if run_opt:
    inp = InputXml()
    inp.read(opt_inp)

    wf_node = inp.find('.//wavefunction')
    j1_node = inp.one_body_jastrow( inp.find_pset(name='ion0') )
    j2_node = inp.two_body_jastrow( inp.find_pset(name='e') )
    for jas_node in [j1_node,j2_node]:
      wf_node.append(jas_node)
    # end for
    inp.write(jas_opt_inp)

    out,err = run_inp_in_mydir(jas_opt_inp,opt_dir,sub_sbatch,pseudo_dir,pseudo_files,name='jopt',submit=True)
    with open(out_dump,'a') as f:
      f.write(out)
      f.write(err)
    # end with
  # end if run_opt

  if run_dmc:

    # get optimized wavefunction
    inp = InputXml()
    inp.read(opt_wf_fname)
    opt_wf_node = inp.find('.//wavefunction')
    # modify h5 ref
    sb_node = opt_wf_node.find('.//sposet_builder')
    sb_node.set('href', os.path.basename(sb_node.get('href')) )

    myinp = InputXml()
    myinp.read(dmc_inp)

    # swap out wavefunction node
    qsys = myinp.find('.//qmcsystem')
    wf_node = qsys.find('.//wavefunction')

    qsys.insert(qsys.index(wf_node),opt_wf_node)
    qsys.remove(wf_node)

    myinp.write(jas_dmc_inp)

    out,err = run_inp_in_mydir(jas_dmc_inp,dmc_dir,sub_sbatch,pseudo_dir,pseudo_files,name='dmc',submit=True)
    with open(out_dump,'a') as f:
      f.write(out)
      f.write(err)
    # end with
  # end if run_opt

# end __main__
