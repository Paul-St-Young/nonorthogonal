#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np
from mmap import mmap

def grab_dets(det_fname,ndet):
  with open(det_fname,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  idx = mm.find('/')
  header = mm[:idx+1]
  header = header.replace('NCI=         800','NCI=         %d'%ndet) # !!!! hard code NCI

  det_idx_list = []
  for idet in range(ndet):
    idx = mm.find('Determinant') 
    if idx == -1:
      raise RuntimeError('%d dets in %s, not enought for requested %d' % (idet,det_fname,ndet))
    # end if
    det_idx_list.append(idx)
    mm.seek(idx)
    mm.readline()
  # end for idet

  idx = mm.find('Determinant') # OK even if idx = -1
  det_text = mm[det_idx_list[0]:idx]

  # fake CI coefficients
  coeff_text = ''
  for idet in range(ndet):
    coeff_text += '(1.0,1.0)\n'
  # end for idet

  return header + "\n" + coeff_text + det_text
# end def grab_dets

if __name__ == '__main__':

  det_fname = '../b-phf/cont_dets/determinants1.det'
  ndet = 25 # desired number of determinants

  get_dir = 'get%d' % ndet
  if not os.path.isdir(get_dir):
    sp.check_call(['mkdir',get_dir])
  # end if

  new_text = grab_dets(det_fname,ndet)
  with open( os.path.join(get_dir,'determinants1.det'),'w' ) as f:
    f.write(new_text)
  # end with

# end __main__
