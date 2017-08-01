#!/usr/bin/env python
import os
import numpy as np

if __name__ == '__main__':
  sub_fname= 'py.sbatch'
  exx2scan = [0.00,0.25,0.50,0.75,1.00]

  for exx in exx2scan:
    subdir = 'exx%d' % (exx*100)

    # copy reference files
    if not os.path.isdir(subdir):
      os.system('cp -r ref %s'%subdir)
    else: # update main.py
      os.system('cp ref/main.py %s/'%subdir)
      os.system('cp ref/launch.py %s/'%subdir)
    # end if

    # edit submission script
    sub_loc = os.path.join(subdir,sub_fname)
    os.system( 'sed -i "s/myname/%s/" %s' % (subdir,sub_loc) )
    os.system( 'sed -i "s/myexx/%3.2f/" %s' % (exx,sub_loc) )
  # end for exx
# end __main__
