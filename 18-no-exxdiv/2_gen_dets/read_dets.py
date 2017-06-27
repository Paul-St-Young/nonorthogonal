#!/usr/bin/env python
import numpy as np
from one_det import read_header,read_next_det
from mmap import mmap

def get_multi_determinants(det_fname,conj_ci=False,conj_det=False):
    with open(det_fname,'r+') as f:
        mm = mmap(f.fileno(),0)
    # end with

    ci_coeff = read_header(mm)
    ndet = len(ci_coeff)

    dets = []
    for idet in range(len(ci_coeff)):
        entry   = read_next_det(mm)
        if entry == -1:
          raise RuntimeError('cannot find the %d out of %d determinant'%(idet,ndet))
        # end if
        idx,det = entry
        if conj_det:
          det = det.T.conj()
        # end if
        dets.append(det)
    # end for idet
    detlist = np.array(dets)

    if conj_ci:
      ci_coeff = ci_coeff.conj()
    # end if

    return ci_coeff,detlist
# end def get_multi_determinants

def save_det(idet):
  fname = 'dets.%d' % idet

  import os
  det_dir = 'coeff%d' % idet
  if not os.path.isdir(det_dir):
    os.system('mkdir %s'%det_dir)
  # end if

  ci_name_fmt  = os.path.join(det_dir,'ci_coeff_conj{conj:d}.dat')
  det_name_fmt = os.path.join(det_dir,'detlist_conj{conj:d}.dat')
  for conj_ci in [False,True]:
    for conj_det in [False,True]:
      ci_coeff,detlist = get_multi_determinants(fname,conj_ci,conj_det)
      assert len(ci_coeff) == len(detlist)
      ci_fname = ci_name_fmt.format(conj=conj_ci)
      det_fname= det_name_fmt.format(conj=conj_det)
      np.savetxt(ci_fname,ci_coeff.view(float))
      np.savetxt(det_fname,detlist.view(float))
    # end for
  # end for
# end def save_det

if __name__ == '__main__':

  idet = 20
  save_det(idet)

# end __main__
"""
      if conj_ci:
        np.savetxt( ci_name_fmt.format(conj=conj_ci), ci_coeff.conj().view(float) )
      else:
        np.savetxt( ci_name_fmt.format(conj=conj_ci), ci_coeff.view(float) )
      # end if
      if conj_ci:
        np.savetxt( det_name_fmt.format(conj=conj_det), detlist.transpose([0,2,1]).conj().view(float) )
      else:
        np.savetxt( det_name_fmt.format(conj=conj_det), detlist.view(float) )
      # end if
""";
