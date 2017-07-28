#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np

def prepare_detdir(subdir,idet,h5_loc,qmc_inp,pseudos):
  if not os.path.isdir(subdir):
    sp.check_call(['mkdir',subdir])
  # end if
  det_dir = os.path.join(subdir,'det%d'%idet) # run location
  if not os.path.isdir(det_dir):
    sp.check_call(['mkdir',det_dir])
  # end if
  inp_loc  = os.path.join(det_dir,qmc_inp)
  h5_href  = os.path.relpath(h5_loc, os.path.dirname(inp_loc) )

  # send pseudopotentials
  for psp in pseudos.itervalues():
    psp_src = os.path.join(pseudo_dir,psp)
    cmd = ' '.join(['cp',psp_src,det_dir])
    #sp.check_call(cmd)
    os.system(cmd)
  # end for
  return det_dir, inp_loc, h5_href
# end def

def jopt_inp(idet,h5_loc,opt_inp,pseudos):
  det_dir, inp_loc, h5_href = prepare_detdir('jopt',idet,h5_loc,opt_inp,pseudos) # setup jopt/det10 for example

  # get <wavefunction> - Jastrow needed!
  coeff_fname = coeff_fname_fmt.format(ndet=idet)
  ci_coeff = np.loadtxt(coeff_fname).view(complex)
  wf_node = inp.uhf_multidet_qmc(ci_coeff,nup,ndn,fftgrid,h5_href=h5_href,real_coeff=False)

  # write QMCPACK input
  inp.write_qmcpack_input(inp_loc,cell,h5_loc,nup,ndn,wf_node=wf_node,pseudos=pseudos,qmc_nodes=[opt_node],proj_id=proj_id)

  # edit QMCPACK input
  ion_cusp = 0.0
  inp.read(inp_loc)
  j1_node = inp.one_body_jastrow( inp.find_pset(name='ion0'),cusp=str(ion_cusp) )
  j2_node = inp.two_body_jastrow( inp.find_pset(name='e') )
  for jas_node in [j1_node,j2_node]:
    wf_node.append(jas_node)
  # end for

  # re-write QMCPACK input
  inp.write_qmcpack_input(inp_loc,cell,h5_loc,nup,ndn,wf_node=wf_node,pseudos=pseudos,qmc_nodes=[opt_node],proj_id=proj_id)
# end def

if __name__ == '__main__':

  idet = 10 # number of determinants in wavefunction

  proj_id = 'mno'
  pseudo_dir = '~/pseudo'
  pseudos = {'Mn':'Mn.BFD.xml','O':'O.BFD.xml'}
  fftgrid = np.array([100,100,100])
  h5_loc  = '../c-dets/pyscf2qmcpack.h5'
  vmc_inp = 'vmc.xml'
  opt_inp = 'opt.xml'
  dmc_inp = 'dmc.xml'
  coeff_fname_fmt = '../d-afqmc/afcoeff{ndet:d}.dat'

  # get simulation cell (nup & ndn)
  chkfile = '../a-hf/bfd.h5'
  from pyscf.pbc import scf
  cell,scf_rec = scf.chkfile.load_scf(chkfile)
  nup,ndn  = cell.nelec

  # get <qmc> section
  from input_xml import InputXml
  inp = InputXml()
  nvmc_walkers = 16   # 16 nodes * 16 walker/node = 256 walkers for OpenMP, 8192 walkers for pure MPI
  ndmc_walkers = 4608 # 16 nodes * 32 cores/node * 8 walkers/core
  nloop = 8
  vmc_node  = inp.get_qmc_node(nvmc_walkers,checkpoint=-1)
  opt_node  = inp.get_optimization_node(nloop)
  dmc_inputs = {
   'substeps':1,
   'steps':10,
   'timestep':1.0,
   'warmupsteps':10,
   'blocks':40,
   'usedrift':'yes'
  }
  dmc_nodes = inp.get_dmc_nodes(ndmc_walkers,nvmc_walkers=nvmc_walkers
   ,correlation_time=0.4,time_step_list=[0.008,0.004],param_name_val_map=dmc_inputs)

  # ============================== run QMCPACK  ============================== 

  # OPT
  #jopt_inp(idet,h5_loc,opt_inp,pseudos)

  # DMC
  jopt_dir = 'jopt'
  opt_jas_fname = 'mno.s007.opt.xml'

  subdir = 'dmc'
  det_dir, inp_loc, h5_href = prepare_detdir(subdir,idet,h5_loc,dmc_inp,pseudos) # setup jopt/det10 for example

  # get <wavefunction> - Jastrow needed!
  coeff_fname = coeff_fname_fmt.format(ndet=idet)
  ci_coeff = np.loadtxt(coeff_fname).view(complex)
  wf_node = inp.uhf_multidet_qmc(ci_coeff,nup,ndn,fftgrid,h5_href=h5_href,real_coeff=False)
  # get <jastrow> 
  jas_loc = os.path.join(jopt_dir,'det%d'%idet,opt_jas_fname)
  jinp = InputXml()
  jinp.read(jas_loc)
  jas_nodes = jinp.findall('.//jastrow')
  for jas in jas_nodes:
    wf_node.append(jas)
  # end for

  # write QMCPACK input
  inp.write_qmcpack_input(inp_loc,cell,h5_loc,nup,ndn,wf_node=wf_node,pseudos=pseudos,qmc_nodes=dmc_nodes,proj_id=proj_id)

# end __main__
