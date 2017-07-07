#!/usr/bin/env python
import numpy as np
def str2comp(text):
  tokens = text.strip(' ()\n').split(',')
  real,imag = map(float,tokens)
  return real+1j*imag
# end def str2comp

def parse_detij(det_fname,idet,jdet,ndet):

  from mmap import mmap
  with open(det_fname,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  # check number of determinants
  idx = mm.find('det =')
  mm.seek(idx)
  det_line = mm.readline() # e.g. MO coefficients [det =   5] (real) :
  deteq = det_line.split('[')[-1].split(']')[0]
  myndet = int(deteq.split('=')[-1])
  assert myndet == ndet

  # check ci coefficients
  idx = mm.find('expansion coefficients')
  mm.seek(idx)
  mm.readline()
  # get the second set of coefficients
  idx = mm.find('expansion coefficients')
  mm.seek(idx)
  mm.readline()
  for i in range(ndet):
    ci_line = mm.readline()
    ci = str2comp(ci_line)
    assert np.isclose(ci.imag,0.0)
    if i==idet or i==jdet:
      assert np.isclose(ci.real,1.0)
    else:
      assert np.isclose(ci.real,0.0)
    # end if
  # end for

  # read energy
  idx = mm.find('hmstat eigenvalues')
  mm.seek(idx)
  # example:
  # hmstat eigenvalues :
  #                        1
  #        1   -7.241590E+00
  #
  mm.readline()
  mm.readline()
  eline = mm.readline()
  tokens = eline.split()
  assert len(tokens) == 2 # make sure only one state is used
  energy = float( eline.split()[1] )
  return energy
# end def parse_detij
  
if __name__ == '__main__':
  ndet = 5
  ewald_correction = -2.695783

  fp = open('paire.dat','w')

  data = []
  from itertools import combinations
  for (idet,jdet) in combinations(range(ndet),2):
    det_fname = 'det%d%d.out' % (idet,jdet)
    edet = parse_detij(det_fname,idet,jdet,ndet)
    entry = [idet,jdet,edet+ewald_correction]
    data.append(entry)
    fp.write('%d  %d %10.6f\n' % tuple(entry) )
  # end for
  fp.close()


# end __main__
