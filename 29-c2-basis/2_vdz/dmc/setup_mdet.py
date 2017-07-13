#!/usr/bin/env python
import os
import numpy as np
from lxml import etree

def setup_folder(dir_name,main_fname,wf_fname,wf_h5_fname,pseudo_file,inp_fname,jas_file=None):

  # get <qmcsystem> from main input template
  parser = etree.XMLParser(remove_blank_text=True)
  doc = etree.parse(main_fname,parser)
  root = doc.getroot()
  qmcsystem = root.find('.//qmcsystem')
  ham_idx = qmcsystem.index( qmcsystem.find('.//hamiltonian') )


  # get <wavefunction> from wavefunction file
  wf = etree.parse(wf_fname,parser).getroot()

  # get <sposet_builder> and change h5 link
  spob = wf.find('.//sposet_builder')
  spob.set('href',wf_h5_fname)

  # complete input
  qmcsystem.insert(ham_idx,wf)

  # append Jastrow to wavefunction if provided
  if jas_file is not None:
    jas = etree.parse(jas_file,parser).getroot()
    jas_list = jas.findall('.//jastrow')
    for jnode in jas_list:
      wf.append(jnode)
    # end for
  # end if

  # write input to local directory
  outfile = os.path.join(dir_name,inp_fname)
  doc.write(outfile,pretty_print=True) 

# end def

def setup_ci_mdet(dir_name,ci_coeff,wf_h5_fname,nfill
  ,pseudo_file = '~/soft/qmcpack/tests/solids/diamondC_1x1x1_pp/C.BFD.xml'
  ,qsub_file = './qmc.msub'
  ,main_fname = '../msd.xml'
  ,jas_file='opt49/c2.s004.opt.xml'
  ,run=False):
  import sys
  #sys.path.insert(0,'../../../19-phfrun-out/3_eigsys')
  from write_mdet import save_mdet
  # make folder if not found
  if not os.path.isdir(dir_name):
    os.system('mkdir %s' % dir_name)
  # end if

  # write wf xml
  wf_fname = os.path.join(dir_name,'mdet.xml')
  save_mdet(wf_fname,ci_coeff,wf_h5_fname,nfill)

  # send pseudopotential
  os.system( 'cp %s %s' % (pseudo_file,dir_name) )

  setup_folder(dir_name,main_fname,wf_fname,wf_h5_fname,pseudo_file,'vmc.xml',jas_file=jas_file)

  # send submission script
  os.system( 'cp %s %s' % (qsub_file,dir_name) )

  subfile = os.path.join(dir_name,os.path.basename(qsub_file))
  os.system('sed -i "s/my_epic_name/%s/" %s' % (dir_name,subfile))

  # !!!! hack submission in here
  if run:
    os.system( 'cd %s; sbatch %s' % (dir_name,os.path.basename(qsub_file)) ) 
# end def

def main():
  import sys
  sys.path.insert(0,'../../../15-pyscf2qmcpack/3_qmcpack')
  from show_fftbox import read_comp
  
  det_name = 'dets'
  wf_h5_fname = '../%s.h5' % det_name
  wf_h5_fname = os.path.join('..',wf_h5_fname) # !!!! hard-code one folder in

  conj_ci = 0

  ndet = 5  # !!!! hard-code
  nfill= 4  # !!!! hard-code
  #ci_coeff = np.loadtxt('../coeff_list/ci%d.dat' % (ndet-1)).view(complex)
  ci_coeff =  np.loadtxt( '../afqmc_ci/afcoeff%d.dat'%(ndet-1) ).view(complex)
  dir_name = det_name+'ci%d'%(ndet-1) + '_af'
  if conj_ci:
    dir_name+='conj'
    ci_coeff = ci_coeff.conj()
  # end if
  setup_ci_mdet(dir_name,ci_coeff,wf_h5_fname,nfill,run=True)

# end def main
  
if __name__ == '__main__':
  main()
# end __main__
