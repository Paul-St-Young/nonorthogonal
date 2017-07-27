#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np

if __name__ == '__main__':

  proj_id = 'mno'
  pseudos = {'Mn':'Mn.BFD.xml','O':'O.BFD.xml'}
  fftgrid = np.array([100,100,100])
  h5_loc  = '../c-dets/pyscf2qmcpack.h5'
  vmc_inp = 'vmc.xml'
  coeff_fname_fmt = '../d-afqmc/afcoeff{ndet:d}.dat'

  # get simulation cell (nup & ndn)
  chkfile = '../a-hf/bfd.h5'
  from pyscf.pbc import scf
  cell,scf_rec = scf.chkfile.load_scf(chkfile)
  nup,ndn  = cell.nelec

  # UHF multi-determinant - need to add Jastrow!
  from input_xml import InputXml
  inp = InputXml()
  nwalker  = 144 # need a lot of walkers to fight variance from cusps
  vmc_node = inp.get_qmc_node(nwalker,checkpoint=-1)

  #   write QMCPACK inputs
  subdir = 'mdet'
  if not os.path.isdir(subdir):
    sp.check_call(['mkdir',subdir])
  # end if

  # ============================== run QMCPACK  ============================== 
  idet = 10 # number of determinants in wavefunction
  det_dir = os.path.join(subdir,'det%d'%idet) # run location
  if not os.path.isdir(det_dir):
    sp.check_call(['mkdir',det_dir])
  # end if
  inp_loc  = os.path.join(det_dir,vmc_inp)
  h5_href  = os.path.relpath(h5_loc, os.path.dirname(inp_loc) )

  # get CI coefficients
  coeff_fname = coeff_fname_fmt.format(ndet=idet)
  ci_coeff = np.loadtxt(coeff_fname).view(complex)

  # write <multideterminant> node
  wf_node  = inp.uhf_multidet_qmc(ci_coeff,nup,ndn,fftgrid,h5_href=h5_href,real_coeff=False)
  
  # write QMCPACK input
  inp.write_qmcpack_input(inp_loc,cell,h5_loc,nup,ndn,wf_node=wf_node,pseudos=pseudos,qmc_nodes=[vmc_node],proj_id=proj_id)

# end __main__
