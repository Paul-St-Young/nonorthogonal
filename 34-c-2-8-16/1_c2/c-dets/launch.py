#!/usr/bin/env python
# launch qmc calculations using template input files written by main.py
import os
import subprocess as sp
from lxml import etree

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

def run_inp_in_mydir(inp_src,mydir,sub_sbatch,pseudo_dir,pseudo_files,submit=True,name='qmc'):
  # make dir
  if not os.path.isdir(mydir):
    sp.check_call(['mkdir',mydir])
  # end if

  # send input
  inp = os.path.basename(inp_src)
  inp_loc = os.path.join(mydir,inp)
  if not os.path.isfile(inp_loc):
    sp.check_call(['cp',inp_src,inp_loc])
  # end if
  doc = etree.parse(inp_loc)
  root = doc.getroot()
  h5_node = root.find('.//sposet_builder')
  wf_h5_fname = h5_node.get('href')
  h5_node.set('href',os.path.join('../../',wf_h5_fname))
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
    queue= 'pdebug',
    minutes = 20,
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

  out = ''; err = ''
  if submit:
    cmd  = 'cd {subdir:s};sbatch {sub_file:s}'.format(subdir=mydir,sub_file=sub_sbatch)
    proc = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
    out,err = proc.communicate()
  # end if submit

  return out,err
# end run_inp_in_mydir

if __name__ == '__main__':

  vmc_inp      = 'vmc.xml'
  pseudo_dir   = '~/pseudo'
  pseudo_files = ['C.BFD.xml']
  sub_sbatch   = 'qmc.sbatch' # name of submision file, written from scratch

  out_dump = 'loc_launch.txt'

  # run everything
  #proc = sp.Popen(' '.join(['ls','-d','./hii/*']),stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
  #out,err = proc.communicate()
  #for subdir in out.split('\n')[:-1]:

  ndet_to_run = 5
  for idet in range(ndet_to_run):
    subdir = './hii/det%d'%idet
    
    job_name = os.path.basename(subdir)
    idet = int( job_name.replace('det','') )

    vmc_src = os.path.join(subdir,vmc_inp)
    out,err = run_inp_in_mydir(vmc_src,subdir,sub_sbatch,pseudo_dir,pseudo_files,name=job_name,submit=True)
    with open(out_dump,'a') as f:
      f.write(out)
      f.write(err)
    # end with
  # end for

# end __main__
