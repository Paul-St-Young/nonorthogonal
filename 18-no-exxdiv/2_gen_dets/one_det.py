#!/usr/bin/env python
import numpy as np
from mmap import mmap

def str2comp(text):
  tokens = text.strip(' ()\n').split(',')
  real,imag = map(float,tokens)
  return real+1j*imag
# end def str2comp

def val(line,sep='=',pos=-1):
    return line.split(sep)[pos].strip()
def get_val(name,mm,sep='=',pos=-1):
    idx = mm.find(name)
    if idx == -1:
        raise RuntimeError('%s not found'%name)
    # end if
    mm.seek(idx)
    line = mm.readline()
    return val(line,sep,pos)

def read_header(mm,header_loc=0):
    mm.seek(header_loc) # default is to rewind file

    # make sure this file is what I expect
    assert mm.find('FullMO') != -1
    assert mm.find('CMajor') != -1
    orb_type = get_val('TYPE',mm)
    assert orb_type == 'rotated'
    mm.seek(header_loc)
    det_format = get_val('FORMAT',mm)
    assert det_format == 'PHF'
    mm.seek(header_loc)
    uhf = int(get_val('UHF',mm))

    # get number of determinants
    mm.seek(header_loc)
    nci = int(get_val('NCI',mm))

    # get coefficients
    idx = mm.find('/')
    mm.seek(idx)
    mm.readline() # skip the line having
    ci_coeff = []
    for ici in range(nci):
      line = mm.readline().strip()
      if line.startswith('(') and line.endswith(')'):
        pair = line.split(',')
        if len(pair) == 2:
          ci_coeff.append(str2comp(line))
        elif len(pair) == 3:
          for num in line.split():
            ci_coeff.append(str2comp(num))
          # end for
        else:
          raise RuntimeError('cannot handle %d numbers on a line'%len(pair))
        # end if
      # end if
    # end for
    assert len(ci_coeff) == nci

    next_line = mm.readline().strip()
    assert next_line.startswith('Determinant:')
    mm.seek(header_loc)

    return np.array(ci_coeff)
# end def


def read_next_det(mm,prefix='Determinant',max_nao=1024):
    """ mm: mmap object
        max_nao: maximum number of atominc orbital in the determinant """

    idx = mm.find(prefix)
    if idx == -1:
      return -1
    mm.seek(idx)

    det_data = []
    det_line = mm.readline()
    det_idx = int( det_line.split(':')[-1] )
    for i in range(max_nao*max_nao):
      cur_idx = mm.tell()
      line = mm.readline().strip()
      if line.startswith('(') and line.endswith(')'):
        pair = line.split(',')
        if len(pair) == 2:
          det_data.append(str2comp(line))
        elif len(pair) == 3:
          for num in line.split():
            det_data.append(str2comp(num))
          # end for
        else:
          raise RuntimeError('cannot handle %d numbers on a line'%len(pair))
        # end if
      else:
        mm.seek(cur_idx)
        break
      # end if
    # end for i

    return det_idx, np.array(det_data,dtype=complex)
# end def read_next_det

if __name__ == '__main__':

    nbas = 48
    det_fname = 'dets.1'

    with open(det_fname,'r+') as f:
      mm = mmap(f.fileno(),0)
    # end with

    idx,data = read_next_det(mm)
    nmo = int(round(np.sqrt(len(data))))
    arrC = data.reshape(nmo,nmo,order='C')
    arrF = data.reshape(nmo,nmo,order='F')
    # arrC is row major, arrF is column major
    assert np.allclose(arrC,arrF.T)

    for icol in range(3,4):
        istart = icol*nbas
        col  = data[istart:istart+nbas]
        #assert np.allclose(col,arrC[icol])
        # consecutive numbers in data belong to the same column
        assert np.allclose(col,arrF[:,icol])
    # end for icol
# end __main__
