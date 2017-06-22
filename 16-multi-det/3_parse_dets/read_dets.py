#!/usr/bin/env python
import numpy as np
from one_det import read_header,read_next_det
from mmap import mmap

def get_multi_determinants(det_fname):
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
        dets.append(det)
    # end for idet
    detlist = np.array(dets)

    return ci_coeff,detlist
# end def get_multi_determinants

if __name__ == '__main__':
    det_fname = '../2_gen_dets/dets.1'
    ci_coeff,detlist = get_multi_determinants(det_fname)
    assert len(ci_coeff) == len(detlist)
    np.savetxt('ci_coeff.dat',ci_coeff.view(float))
    np.savetxt('detlist.dat',detlist.view(float))
# end __main__
