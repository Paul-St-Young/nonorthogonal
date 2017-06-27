#!/usr/bin/env python
import os
import numpy as np
from lxml import etree

def setup_folder(dir_name,main_fname,wf_fname,wf_h5_fname,pseudo_file,inp_fname):

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

  # complete input and write to local directory
  qmcsystem.insert(ham_idx,wf)
  outfile = os.path.join(dir_name,inp_fname)
  doc.write(outfile,pretty_print=True) 

# end def

def main():
  import sys
  sys.path.insert(0,'../3_eigsys')
  from write_mdet import save_mdet

  main_fname = '../4_run/msd.xml'
  pseudo_file = '~/soft/qmcpack/pseudopotentials/BFD/C.BFD.xml'
  qsub_file = '~/notes/qmc.msub'
  wf_h5_fname = '../3_eigsys/dets.h5'
  wf_h5_fname = os.path.join('..',wf_h5_fname) # !!!! hard-code one folder in

  for idet in range(0,8):
    ndet = 21 # !!!! hard-code
    nfill= 4  # !!!! hard-code
    ci_coeff = np.zeros(ndet)
    ci_coeff[idet] = 1.0

    # write wf xml
    wf_fname = 'mdet%d.xml' % idet
    save_mdet(wf_fname,ci_coeff,wf_h5_fname,nfill)
    dir_name = 'det%d'%idet

    # make folder if not found
    if not os.path.isdir(dir_name):
      os.system('mkdir %s' % dir_name)
    # end if
    # send pseudopotential
    os.system( 'cp %s %s' % (pseudo_file,dir_name) )
    setup_folder(dir_name,main_fname,wf_fname,wf_h5_fname,pseudo_file,'vmc.xml')

    # send submission script
    os.system( 'cp %s %s' % (qsub_file,dir_name) )

    subfile = os.path.join(dir_name,os.path.basename(qsub_file))
    os.system('sed -i "s/my_epic_name/%s/" %s' % (dir_name,subfile))

    # !!!! hack submission in here
    os.system( 'cd %s; sbatch %s' % (dir_name,os.path.basename(qsub_file)) ) 

    # end for
  # end for
# end def main
  
if __name__ == '__main__':
  main()
# end __main__
