#!/usr/bin/env python
# launch qmc calculations using template input files written by main.py
import os
import subprocess as sp
from lxml import etree
from input_xml import InputXml

def submission_template():
  ftext = '''#!/bin/bash
#SBATCH --job-name {name:s}
#SBATCH -N {nodes:d}
#SBATCH -p {queue:s}
#SBATCH -t {minutes:d}
#SBATCH --account qmchhp

BIN={bin_loc:s}

export OMP_NUM_THREADS={threads:d}
date
srun -n {nproc:d} $BIN {inp:s} > out 2> err
date '''
  return ftext
# end def submission_template

def run_inp_in_mydir(inp,mydir,sub_sbatch,pseudo_dir,pseudo_files,submit=True,name='qmc'):
  nodes = 4
  minutes= 120
  queue = 'pbatch'
  account = 'qmchhp'
  npp   = 36
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
  nproc = nodes*npp
  sub_content = sub_text.format(
    name = name,
    account = account,
    nodes= nodes,
    nproc= nproc,
    queue= queue,
    minutes = minutes,
    bin_loc = '/g/g91/yang41/soft/master_qmcpack/build/bin/qmcpack',
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

  run_vmc = False
  run_opt = False
  run_dmc = True
  opt_wf_fname = 'jopt/mno.s004.opt.xml'

  pseudo_dir   = '~/pseudo'
  pseudo_files = ['Mn.BFD.xml','O.BFD.xml']
  ion_cusp = 0.0 # !!!! do not add cusp to electron-ion Jastrow with pseudo pot.
  sub_sbatch   = 'qmc.sbatch' # name of submision file, written from scratch

  vmc_dir = 'vmc'
  vmc_inp = 'vmc.xml'

  opt_dir = 'jopt'
  opt_inp = 'opt.xml'
  jas_opt_inp = 'jopt.xml'

  dmc_dir = 'dmc'
  dmc_inp = 'dmc.xml'
  jas_dmc_inp = 'jdmc.xml'

  out_dump = 'loc_launch.txt'

  if run_vmc:
    out,err = run_inp_in_mydir(vmc_inp,vmc_dir,sub_sbatch,pseudo_dir,pseudo_files,name='vmc',submit=True)
    with open(out_dump,'a') as f:
      f.write(out)
      f.write(err)
    # end with
  # end if

  if run_opt:
    inp = InputXml()
    inp.read(opt_inp)

    wf_node = inp.find('.//wavefunction')
    j1_node = inp.one_body_jastrow( inp.find_pset(name='ion0'),cusp=str(ion_cusp) )
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
