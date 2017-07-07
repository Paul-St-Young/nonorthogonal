#!/usr/bin/env python
import os
import numpy as np
from lxml import etree

def setup_folder(dir_name,main_fname,wf_h5_fname,inp_fname):

  # get <qmcsystem> from main input template
  parser = etree.XMLParser(remove_blank_text=True)

  # get <wavefunction> from wavefunction file
  doc = etree.parse(main_fname,parser)
  wf  = doc.getroot()

  # get <sposet_builder> and change h5 link
  spob = wf.find('.//sposet_builder')
  spob.set('href',wf_h5_fname)

  # complete input and write to local directory
  outfile = os.path.join(dir_name,inp_fname)
  doc.write(outfile,pretty_print=True) 

# end def

def main():

  main_fname = './ref.xml'
  qsub_file = './qmc.msub'

  ndet = 5 # !!!! hard-code
  for idet in range(0,ndet):
    wf_h5_fname = '../det%d.h5' % idet
    wf_h5_fname = os.path.join('..',wf_h5_fname) # !!!! hard-code one folder in

    dir_name = 'det%d'%idet

    # make folder if not found
    if not os.path.isdir(dir_name):
      os.system('mkdir %s' % dir_name)
    # end if

    setup_folder(dir_name,main_fname,wf_h5_fname,'vmc.xml')

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
